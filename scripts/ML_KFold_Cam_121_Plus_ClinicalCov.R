#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

classifierName = args[1]
group = args[2]
pathForTrainDataSet = args[3]
pathForTestDataSet = args[4]
pathForGeneToGroupMappings = args[5]
applyPCA = args[6]
pathForHelperFile = args[7]
outputFileNameForFinalResults = args[8]
outputFileNameForTrainPredBestTune = args[9]
outputFileNameForTestPredBestTune = args[10]
seed = args[11]
pathForNumPcs = args[12]

# Load the Gene to Group mappings
geneToGroupMappings = read.table(pathForGeneToGroupMappings, sep = "\t", header = TRUE, quote = "")
source(pathForHelperFile)
extractedData = extractExpressionAndClinicalData(pathForTrainDataSet, geneToGroupMappings, group)

rawCountsExpressionData = extractedData[[1]]
clinicalData = extractedData[[2]] 

############################################ TRAIN ############################################
#Pre-process the data
## Normalize the expression data to stablize the variance
#Apply variance stabilizing transformation on the training data as explained here: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation
#colData = apply(clinicalData, 2, as.factor)
ddsTrain = DESeqDataSetFromMatrix(countData = t(rawCountsExpressionData), colData = clinicalData, design = ~Stage+Bres+ECOG+treatment+EventMet)
# learn the dispersion function of a dataset
ddsTrain <- estimateSizeFactors(ddsTrain)
ddsTrain <- estimateDispersions(ddsTrain)
xTrainExprs <- data.frame(t(assay(varianceStabilizingTransformation(ddsTrain, blind = FALSE))))

#Pre-process the clinical data
#Converting every categorical variable to numerical using dummy variables
dmy <- dummyVars(" ~ .", data = clinicalData[, c("Stage", "Bres", "ECOG", "treatment")],fullRank = F) #We would need dummy variables for all the factor levels
xTrainCov <- data.frame(predict(dmy, newdata = clinicalData[, c("Stage", "Bres", "ECOG", "treatment")]))
#Removing the redundant columns
xTrainCov <- subset(xTrainCov, select = -c(ECOG.Asymptomatic,treatment.Observation))

yTrain <- clinicalData$EventMet


#Find the number of principal components to retain if applyPCA = TRUE
if(applyPCA){
  #Automatically determine the number of PCs using the ElbpwFinder function from the package "RclusTool"
  #library("RclusTool")
  #pcaResults = prcomp(xTrain)
  #numPCs = as.numeric(ElbowFinder(y = pcaResults$sdev^2, x = 1:length(pcaResults$sdev)))
  #df = read.table(pathForNumPcs, sep = "\t", header = TRUE, quote = "")
  #numPCs = as.numeric(df$numPCs[df$group == group])
  numPCs = 17
} else{
  numPCs = 0
}

#Train the model
##Combine the expression data and the clinical data
xTrain = cbind(xTrainExprs, xTrainCov)
resultsWithTrainedModels = trainModel(entireData = xTrain, entireLabels = yTrain, seed = seed, applyPCA = applyPCA, pcaComp = numPCs, classifierName = classifierName)

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
############################################ TEST ############################################
## Load the testing data
## Extract the raw counts
extractedData = extractExpressionAndClinicalData(pathForTestDataSet, geneToGroupMappings, group)

rawCountsExpressionData = extractedData[[1]]
clinicalData = extractedData[[2]] 
## Make sure that the variable names are in the same order as the training data
rawCountsExpressionData = rawCountsExpressionData[, colnames(xTrainExprs)]

## Preprocess the expression data
ddsTest = DESeqDataSetFromMatrix(countData = t(rawCountsExpressionData), colData = clinicalData[!colnames(clinicalData)%in%c("EventMet")], design = ~1) #Exclude the "EventMet" column
ddsTest = estimateSizeFactors(ddsTest)
dispersionFunction(ddsTest) = dispersionFunction(ddsTrain)
xTestExprs = data.frame(t(assay(varianceStabilizingTransformation(ddsTest, blind = FALSE))))

#Preprocess the clinical data
xTestCov <- data.frame(predict(dmy, newdata = clinicalData[, c("Stage", "Bres", "ECOG", "treatment")]))
#Remove the redundant columns
xTestCov <- subset(xTestCov, select = -c(ECOG.Asymptomatic,treatment.Observation))

yTest<-clinicalData$EventMet

## Predict performance of the final model on the test dataset
###Combine the expression data and the clinical data
xTest = cbind(xTestExprs, xTestCov)
predTest = predict(resultsWithTrainedModels, xTest, type = "prob")$model
predTest$obs = yTest
############################################ SAVE ############################################
write.table(predTest, file = outputFileNameForTestPredBestTune, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

