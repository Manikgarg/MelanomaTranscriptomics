library("DESeq2")
#Install apeglm (v: "1.6.0") for R (v "3.6")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("apeglm")
library("apeglm")
fileNames = list.files(path = "./Desktop", pattern=".deseq2", full.names=TRUE)
#load(paste(fileNames[1],"de.Rdata", sep = "/"))
#resultsNames(dds)
#resWithLfcShrinkage = lfcShrink(dds, coef = "EventMet_Yes_vs_No", type = "apeglm")

performLfcShrink = function(path){
  #Load the saved DESeq2 results
  load(paste(path,"de.Rdata", sep = "/"))
  #Apply lfc shrinkage
  coefName = resultsNames(dds)
  resWithLfcShrinkage = lfcShrink(dds, coef = tail(coefName, n=1), type="apeglm")
  #Add suffix to colnames for easy identification later on.
  colnames(resWithLfcShrinkage) = paste(colnames(resWithLfcShrinkage),"lfcShrinkApplied", sep = "_")
  #Add ENSEMBL ID as a column in the data frame to ease merging later
  resWithLfcShrinkage$id = rownames(resWithLfcShrinkage)
  
  #Append the new results to old results.
  #Load the old results
  #deResults = read.table(file = paste(path,"de.tsv", sep = "/"), header = TRUE, sep = '\t')
  #Append new results to old results by ENSEMBL ID. 
  #Note that the old results can also be accessed by res.annot, therefore, no need to read the de.tsv file
  deResultsUpdated = merge(x = res.annot, y = data.frame(resWithLfcShrinkage), by = "id")
  #Order the results in increasing order of padj
  deResultsUpdated = deResultsUpdated[order(deResultsUpdated$padj),]
  #Save the updated results to a new file.
  write.table(deResultsUpdated, file = paste(path,"deWithLfcSkrinkageApplied.tsv", sep = "/"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  
  #Create ranked list for running Pre-ranked GSEA tool
  rankedList = na.omit(data.frame(deResultsUpdated$Name, deResultsUpdated$log2FoldChange_lfcShrinkApplied))
  rankedList = rankedList[with(rankedList, order(deResultsUpdated.log2FoldChange_lfcShrinkApplied)), ]
  #Save the ranked list to a new file for GSEA.
  write.table(rankedList, file = paste(path,"rankedList.rnk", sep = "/"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

}

for (path in fileNames){
  print(path)
  #print(paste(path,"de.Rdata", sep = "/"))
  performLfcShrink(path)
}
