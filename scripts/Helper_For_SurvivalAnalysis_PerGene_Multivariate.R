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
  rownames(data) = make.names(rownames(data))
  
  return(list(data, clinicalData))
}


saveSurvivalResults = function(survivalModel){ 
  x <- summary(survivalModel)
  p.value<- signif(x$coef[1, ncol(x$coef)], digits = 5)
  HR <- signif(x$coef[1, "exp(coef)"], digits = 5)
  HR.confint.lower <- signif(x$conf.int[1,"lower .95"], digits = 5)
  HR.confint.upper <- signif(x$conf.int[1,"upper .95"], digits = 5)
  HR <- paste0(HR, " (", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res<-c(HR, p.value)
  names(res)<-c("HR (95% CI for HR)", 
                "p.value")
  return(res)
}
