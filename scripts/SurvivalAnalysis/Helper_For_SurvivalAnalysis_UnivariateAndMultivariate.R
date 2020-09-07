suppressPackageStartupMessages(library(DESeq2))

extractExpressionAndClinicalDataSet = function(pathToTheDESeq2Results, geneToGroupMappings, group, pathToSurvivalMetadata){
  #Load the DESeq2 results
  load(paste0(pathToTheDESeq2Results))
  #Subset the VST normalized expression data to the genes belong to a particular "group"
  select = rownames(res.annot)[ which(res.annot$Name %in% geneToGroupMappings$Gene[geneToGroupMappings[[group]] == 1])]  
  #data = t(data.frame(data.f[select, ]))
  #Taking transpose of the data because LMC_150_genes have duplicated gene names and will throw an error in subsequent analysis. Hence taking double transpose of the data to bypass that error.
  data = t(data.frame(assay(vsd)[select, ]))
  #Replace ENSEMBL IDs with corresponding gene names
  colnames(data) = res.annot$Name[which(rownames(res.annot) %in% select)]
  
  #Remove genes with duplicated names
  data = t(data[, !duplicated(colnames(data))])
  dim(data)
  
  #Extract the clinical data
  clinicalData = data.frame(colData(vsd))
  
  #Make sure that the clinical data is in the same order as samples in the expression data
  seIdx = match(rownames(clinicalData), colnames(data))
  data = data[, seIdx]
  
  return(list(data, clinicalData))
}


# extractExpressionAndClinicalDataSet = function(pathToTheDESeq2Results, geneToGroupMappings, group, pathToSurvivalMetadata){
#   #Load the DESeq2 results
#   load(paste0(pathToTheDESeq2Results))
#   #Extract the expression data
#   expressionData = assay(vsd)
#   #Convert ENSEMBL IDs to gene names
#   seIdx = match(rownames(res.annot), rownames(expressionData))
#   rownames(expressionData) = make.names(res.annot$Name[seIdx])
#   
#   expressionData = expressionData[make.names(geneToGroupMappings$Gene[geneToGroupMappings[[group]] == 1]), ]
#   
#   #Extract the clinical data
#   clinicalData = data.frame(colData(vsd)[, c("Age", "Stage", "Bres", "ECOG", "treatment", "Nclass", "EventMet")])
#   clinicalData$Age = as.numeric(clinicalData$Age)
#   clinicalData$Bres = as.factor(clinicalData$Bres)
#   clinicalData$ECOG = as.factor(clinicalData$ECOG)
#   clinicalData$treatment = as.factor(clinicalData$treatment)
#   clinicalData$Nclass = as.factor(clinicalData$Nclass)
#   clinicalData$Stage = as.factor(clinicalData$Stage)
#   
#   #Add survival data including distant metastases free survival time and relapse free survival time. Referring to the script sent by Johannes based on the metadata provided by us
#   # load the data
#   data <- read.table(pathToSurvivalMetadata, sep = "\t", header = TRUE)
#   data$time_study_entry <- as.numeric(
#     gsub("\\([ ]*([0-9.]*),[ ]*([0-9]*).*", "\\1", as.character(data$survival_d_ltrc)))
#   data$time_until_death <- as.numeric(
#     gsub("\\([ ]*([0-9.]*),[ ]*([0-9]*).*", "\\2", as.character(data$survival_d_ltrc)))
#   data$time_until_progression <- as.numeric(
#     gsub("\\([ ]*([0-9.]*),[ ]*([0-9]*).*", "\\2", as.character(data$survival_rd_ltrc)))
#   
#   #Add these columns to my clinical data
#   #clinicalData = clinicalData[which(rownames(clinicalData) %in% rownames(data)), ]
#   clinicalData$total_time_until_death <- as.numeric(data$time_study_entry + data$time_until_death)
#   clinicalData$total_time_until_progression <- as.numeric(data$time_study_entry + data$time_until_progression)
#   clinicalData$Dead <- data$Dead
#   clinicalData$PFSEvent <- ifelse((data$Dead == FALSE & clinicalData$EventMet == "No"), FALSE, TRUE)
#   
#   #Remove samples with missing information in any of the covariates
#   clinicalData = clinicalData[complete.cases(clinicalData), ]
#   
#   #Make sure that the clinical data is in the same order as samples in the expression data
#   seIdx = match(rownames(clinicalData), colnames(expressionData))
#   expressionData = expressionData[, seIdx]
#   
#   return(list(data, clinicalData))
# }


applyMethod = function(expressionData, method, pathToTheDESeq2Results, pathForBetaCoeff){
  if (method == "High") {
    reducedSignature = apply(expressionData, 2, median)
  }
  else if (method == "PCA1"){
    pcaResults = prcomp(t(expressionData))
    reducedSignature = pcaResults$x[, 1]
  }
  else if (method == "Weighted"){
    #Load the DESeq2 results
    load(paste0(pathToTheDESeq2Results))
    betaCoeff = read.table(pathForBetaCoeff, sep = "\t", header = TRUE, quote = "")
    #Make sure that the clinical data is in the same order as samples in the expression data
    seIdx = match(rownames(res.annot), rownames(betaCoeff))
    betaCoeff = betaCoeff[seIdx, ]
    betaCoeff$Name = res.annot$Name
    betaCoeffSubset = betaCoeff[which(betaCoeff$Name %in% rownames(expressionData)), ]
    seIdx = match(rownames(expressionData), betaCoeffSubset$Name)
    betaCoeffSubset = betaCoeffSubset[seIdx, ]
    reducedSignature = apply(expressionData*betaCoeffSubset$EventMet_Yes_vs_No, 2, sum)
    stand.fun = function(x){(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)}
    reducedSignature = stand.fun(reducedSignature)  
}
  else{
    print("Method not supported")
  }
  return(reducedSignature)
}

divideIntoCategories = function(reducedSignature, threshold){
  threshold = ifelse(threshold == "median", median(reducedSignature), as.numeric(quantile(reducedSignature,as.numeric(threshold))))
  categories = ifelse(reducedSignature >= threshold, "High", "Low")
  return(categories)
}

saveSurvivalResults = function(survivalModel){ 
  x <- summary(survivalModel)
  p.value<- signif(x$coef[1, ncol(x$coef)], digits = 5)
  HR <- signif(x$coef[1, "exp(coef)"], digits = 5)
  HR.confint.lower <- signif(x$conf.int[1,"lower .95"], digits = 5)
  HR.confint.upper <- signif(x$conf.int[1,"upper .95"], digits = 5)
  res<-c(HR, HR.confint.lower, HR.confint.upper, p.value)
  names(res)<-c("HazardRatio","Lower_HazardRatio","Upper_HazardRatio", "p-value")
  return(res)
}
