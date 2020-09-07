#[Train:Validation split = 4:1 and Train:Test split in CV = 4:1 using 5-fold CV, repeated 5 times]repeated 100 times with different seeds
#Scaling the train and test data before training a classifier to have zero mean and unit variance
#24.07.2019: Include Nclass_binary in the list of covariates
#clinicalData = colData(vsd)[seIdx, c("Tissue_Code", "Stage", "Bres", "ECOG", "treatment", "EventMet")]
#Changes clinical data to data.frame
#Change the clinical factor to supervised factors to be sure that different categories get consistent numbers throughout the analysis:
#For example: Stage = "IIIA" always gets the number 3 and ECOG = "Symptomatic" always gets the number 2
#Alex Jung: Factor() can be used to specify levels and also whether the variable is ordinal or nominal. 
#It can also apply the same factor information to the new dataset from the clinicians as long as they use the same notation

#Use the same standardization obtained from training to the testing instead of using scale function separately on training and testing dataset
#i.e. follow point 3 of the following response: https://stats.stackexchange.com/questions/174823/how-to-apply-standardization-normalization-to-train-and-testset-if-prediction-i

#Fix gene names using make.names()
suppressMessages(library("DESeq2"))
suppressMessages(library("caret"))
suppressMessages(library("doParallel"))

extractExpressionAndClinicalData = function(pathToTheDESeq2Results, geneToGroupMappings, group){
  
  #Load the DESeq2 results
  load(paste0(pathToTheDESeq2Results))
  
  if (group != "ClinicalData") {
    #Subset the VST normalized expression data to the genes belong to a particular "group"
    select = rownames(res.annot)[ which(res.annot$Name %in% geneToGroupMappings$Gene[geneToGroupMappings[[group]] == 1])]
    data = t(data.frame(data.f[select, ]))
    #data = t(data.frame(assay(vsd)[select, ]))
    #Replace ENSEMBL IDs with corresponding gene names
    colnames(data) = res.annot$Name[which(rownames(res.annot) %in% select)]
    
    #Remove genes with duplicated names
    data = data[, !duplicated(colnames(data))]
    dim(data)
    
    #Extract the clinical data
    clinicalData = data.frame(colData(vsd)[, c("Stage", "Bres", "ECOG", "treatment", "EventMet")])
    #Remove samples with missing information in any of the covariates
    clinicalData = clinicalData[complete.cases(clinicalData), ]
    
    #Make sure that the clinical data is in the same order as samples in the expression data
    seIdx = match(rownames(clinicalData), rownames(data))
    data = data[seIdx, ]
    colnames(data) = make.names(colnames(data))
    
    return(list(data, clinicalData))
  } else{
    #Extract the clinical data
    clinicalData = data.frame(colData(vsd)[, c("Stage", "Bres", "ECOG", "treatment", "EventMet")])
    #Remove samples with missing information in any of the covariates
    clinicalData = clinicalData[complete.cases(clinicalData), ]
    return(list(clinicalData))
  }
  
}

trainModel = function(entireData, entireLabels, seed, applyPCA = FALSE, pcaComp = 2, classifierName = "rf"){
  #PreProcess the data
  #Check if PCA is requested
  methodsForPreprocessing = (if (applyPCA) c("pca") else c("nzv") )
  
  set.seed(123)
  seeds <- vector(mode = "list", length = 10001)
  for(i in 1:10000) seeds[[i]] <- sample.int(100000, 100)
  
  ## For the last model:
  seeds[[10001]] <- sample.int(10000, 1)
  
  control <- trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 1000,
                          seeds = seeds,
                          search = "random", 
                          preProcOptions = list(pcaComp = pcaComp),
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary,
                          savePredictions = TRUE)
  
  #setup parallel backend to use many processors
  cores = 10
  cl = makeCluster(cores, type = "SOCK") #or cores -1 : not to overload your computer
  registerDoParallel(cl)
  
  set.seed(seed)
  fit.classifier <- train(x = entireData, y = entireLabels, method= classifierName, trControl=control,
                          tuneLength = 100, 
                          preProcess = methodsForPreprocessing)
  
  ## When you are done:
  stopCluster(cl)
  
  return(list(model = fit.classifier))
}
