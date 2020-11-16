args = commandArgs(trailingOnly=TRUE)

classifierName = args[1]
group = args[2]
pathForTrainDataSet = args[3]
#pathForTestDataSet = args[4]
pathForGeneToGroupMappings = args[4]
applyPCA = args[5]
pathForHelperFile = args[6]
outputFileNameForFinalResults = args[7]
outputFileNameForTrainPredBestTune = args[8]
#outputFileNameForTestPredBestTune = args[9]
seed = args[9]
pathForNumPcs = args[10]
cvMethod = args[11]

# Load the Gene to Group mappings
geneToGroupMappings = read.table(pathForGeneToGroupMappings, sep = "\t", header = TRUE, quote = "")
source(pathForHelperFile)
extractedData = extractExpressionAndClinicalData(pathForTrainDataSet, geneToGroupMappings, group)
if (group!="ClinicalData") {
  rawCountsExpressionData = extractedData[[1]]
  clinicalData = extractedData[[2]] 
}else{
  clinicalData = extractedData[[1]]
}

############################################ TRAIN ############################################
#Pre-process the data
if (group!="ClinicalData") {
  ## Normalize the expression data to stablize the variance
  #Apply variance stabilizing transformation on the training data as explained here: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation
  #colData = apply(clinicalData, 2, as.factor)
  ddsTrain = DESeqDataSetFromMatrix(countData = t(rawCountsExpressionData), colData = clinicalData, design = ~Stage+Bres+ECOG+treatment+EventMet)
  # learn the dispersion function of a dataset
  ddsTrain <- estimateSizeFactors(ddsTrain)
  ddsTrain <- estimateDispersions(ddsTrain)
  xTrain <- data.frame(t(assay(varianceStabilizingTransformation(ddsTrain, blind = FALSE))))
}else{
  #Pre-process the clinical data
  #Converting every categorical variable to numerical using dummy variables
  dmy <- dummyVars(" ~ .", data = clinicalData[, c("Stage", "Bres", "ECOG", "treatment")],fullRank = F) #We would need dummy variables for all the factor levels
  xTrain <- data.frame(predict(dmy, newdata = clinicalData[, c("Stage", "Bres", "ECOG", "treatment")]))
  #Removing the redundant columns
  xTrain <- subset(xTrain, select = -c(ECOG.Asymptomatic,treatment.Observation))
}

yTrain <- clinicalData$EventMet


#Find the number of principal components to retain if applyPCA = TRUE
if(applyPCA){
  #Automatically determine the number of PCs using the ElbpwFinder function from the package "RclusTool"
  #library("RclusTool")
  #pcaResults = prcomp(xTrain)
  #numPCs = as.numeric(ElbowFinder(y = pcaResults$sdev^2, x = 1:length(pcaResults$sdev)))
  df = read.table(pathForNumPcs, sep = "\t", header = TRUE, quote = "")
  numPCs = as.numeric(df$numPCs[df$group == group])

} else{
  numPCs = 0
}

#Train the model
resultsWithTrainedModels = trainModel(entireData = xTrain, entireLabels = yTrain, cvMethod = cvMethod, seed = seed, applyPCA = applyPCA, pcaComp = numPCs, classifierName = classifierName)

#Save the model params for interpreting and comparing the results
bestTuneParams = resultsWithTrainedModels[["model"]][["bestTune"]]
if (bestTuneParams != "none" && ncol(bestTuneParams) > 1) {
  m = resultsWithTrainedModels[["model"]][["pred"]]
  bestTuneParams = resultsWithTrainedModels[["model"]][["bestTune"]]
  indices = apply(m[, names(bestTuneParams)], 1, function(r) any(r %in% bestTuneParams))
  predTrain = m[indices, ]
  
  m1 = resultsWithTrainedModels[["model"]][["results"]]
  indices = apply(m1[, names(bestTuneParams)], 1, function(r) any(r %in% bestTuneParams))
  resultsTrain = m1[indices, ]

} else if (bestTuneParams != "none" && ncol(bestTuneParams) == 1) {
  m = resultsWithTrainedModels[["model"]][["pred"]]
  bestTuneParams = resultsWithTrainedModels[["model"]][["bestTune"]]
  indices = which(m[, names(bestTuneParams)] %in% bestTuneParams, arr.ind = TRUE)
  predTrain = m[indices, ]
  
  m1 = resultsWithTrainedModels[["model"]][["results"]]
  indices = which(m1[, names(bestTuneParams)] %in% bestTuneParams, arr.ind = TRUE)
  resultsTrain = m1[indices, ]
  
} else{
  predTrain = resultsWithTrainedModels[["model"]][["pred"]]
  resultsTrain = resultsWithTrainedModels[["model"]][["results"]]
}

############################################ SAVE ############################################
write.table(resultsTrain, file = outputFileNameForFinalResults, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(predTrain, file = outputFileNameForTrainPredBestTune, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#Save the trained ML model in case of rf and plr as the values may need to be extracted manually
if (classifierName == "plr" || classifierName == "rf"){
  modelName = paste("trainedModel", group, applyPCA, classifierName, ".Rdata", sep = "_")
  save(resultsWithTrainedModels, file = modelName)
}
