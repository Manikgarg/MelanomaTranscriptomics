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
    #seIdx = match(rownames(res.annot), rownames(betaCoeff))
    #betaCoeff = betaCoeff[seIdx, ]
    #betaCoeff$Name = res.annot$Name
    betaCoeffSubset = betaCoeff[which(rownames(betaCoeff) %in% rownames(expressionData)), ]
    seIdx = match(rownames(expressionData), rownames(betaCoeffSubset))
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
  p.value<-signif(x$coef[1, ncol(x$coef)], digits=2)
  HR <-signif(x$coef[1, "exp(coef)"], digits=2);
  HR.confint.lower <- signif(x$conf.int[1,"lower .95"], digits = 2)
  HR.confint.upper <- signif(x$conf.int[1,"upper .95"], digits = 2)
  HR <- paste0(HR, " (", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res<-c(HR, p.value)
  names(res)<-c("HR (95% CI for HR)", 
                "p.value")
  return(res)
}

.sigCheckPval = function(performance,scores,lt=T) {
  if(lt) {
    pval <- 1 - (sum(performance < scores)/length(scores))
  } else {
    pval <- 1 - (sum(performance > scores)/length(scores))
  }
  accuracy <- ceiling(log10(length(scores)))
  pval = signif(pval,accuracy)
  return(pval)
}

plotDensity <- function(survivalPval, randomResult, titlestring, fileName){
  #Run the function for calculating p-values from sigcheck package.
  allScoreValues <- c(survivalPval,randomResult)
  scores     <- sort(allScoreValues, decreasing=FALSE)
  geneListRank   <- mean(which(scores == survivalPval))
  checkPval      <- .sigCheckPval(survivalPval,randomResult,lt=TRUE)
  #accuracy <- ceiling(log10(length(sortScores)))
  #modePerformance = ifelse((metricName == "Accuracy"), 0.538, 0.20)
  
  xmin <- min(scores)
  xmax <- max(scores)  
  pval <- 1 - (sum(survivalPval < scores)/length(scores))
  accuracy <- ceiling(log10(length(scores)))
  
  if(pval > 0) {
    rankstr <- sprintf("Percentile:%%.2f (Tests:%%d  p=%%1.%df)",
                       accuracy)
    pvalstr <- sprintf("Significant p-val =%%1.%df",accuracy)
  } else {
    rankstr <- sprintf("Percentile:%%.2f (Tests:%%d  p<%%1.%df)",
                       accuracy)
    pvalstr <- sprintf("Signature   p-val <%%1.%df",accuracy)        
    pval <- 1/(10^accuracy)
  }
  
  if(survivalPval > 0.001) {
    sigpvalstr <- "Signature   p-val =%1.3f (%1.2f)"
    pvalSurvival <- survivalPval
  } else {
    sigpvalstr <- "Signature   p-val <%1.3f (%1.2f)"        
    pvalSurvival <- 1/(10^3)
  }
  
  rankstr <- sprintf(rankstr,1-ecdf(scores)(survivalPval),
                     length(scores),pval)
  
  pdf(file = fileName)
  plot(density(-log10(scores)),
       main=titlestring,sub=rankstr)
  abline(v=-log10(survivalPval),col="red",lwd=3)
  abline(v=-log10(.05),col="red",lty="dotted",lwd=2)
  legend("topright",legend=c(sprintf(sigpvalstr,pvalSurvival, 
                                     -log10(survivalPval)),
                             "Significant p-val <0.05    (1.30)"),
         col="red",lty=c("solid","dotted"),lwd=c(3,2))
  dev.off()
}
  
