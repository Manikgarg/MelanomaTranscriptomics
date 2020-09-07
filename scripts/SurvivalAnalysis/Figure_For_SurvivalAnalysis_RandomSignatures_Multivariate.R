#!/usr/bin/env Rscript
.libPaths("/homes/manikg/R/x86_64-redhat-linux-gnu-library/3.5")
args = commandArgs(trailingOnly=TRUE)

pathToResults = args[1]
pathToHelperFile = args[2]
valueForOS = args[3]
valueForPFS = args[4]

pattern = "\\d{0}_Weighted_[0-9].*tsv$"
tempFiles = list.files(path = pathToResults, pattern = pattern, recursive = FALSE)
randomResult = data.frame()
for (temp in tempFiles) {
  tempResults = read.table(file = paste0(pathToResults, temp), header = TRUE, quote = "", sep = "\t")
  
  #Include other parameters 
  fileName = unlist(strsplit(temp, split='/', fixed=TRUE))
  tempFileName = unlist(strsplit(fileName[length(fileName)], split='_', fixed=TRUE))
  tempResults$SignatureLength = tempFileName[2]
  tempResults$ScoreThreshold = tempFileName[5]
  tempResults$ScoreMethod = tempFileName[4]
  tempResults$SurvivalType = tempFileName[3]
  tempResults$Seed = unlist(strsplit(tempFileName[length(tempFileName)], split='.', fixed=TRUE))[1]
  randomResult = rbind(randomResult, tempResults)
}

source(pathToHelperFile)
for (k in unique(randomResult$SurvivalType)) {
          originalPval = ifelse(k == "PFS", as.numeric(valueForPFS), as.numeric(valueForOS))
          randomResultSelected = randomResult[(#randomResult$SignatureLength == i & 
                                         #randomResult$Background == j &
                                         randomResult$SurvivalType == k #&
                                         #randomResult$ScoreMethod == l &
                                         #randomResult$ScoreThreshold == m
                                         ), "p.value"]
          
          #metricName = "PVal"
          titlestring = paste(k, ": Random Genes")
          fileName = paste("RandomResult", k, sep = "_")
          plotDensity(originalPval, randomResultSelected, titlestring = titlestring, paste0(pathToResults,fileName, ".pdf"))
        }
