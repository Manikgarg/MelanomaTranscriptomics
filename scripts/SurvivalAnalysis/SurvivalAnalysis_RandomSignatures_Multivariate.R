#!/usr/bin/env Rscript
.libPaths("/homes/manikg/R/x86_64-redhat-linux-gnu-library/3.5")

args = commandArgs(trailingOnly=TRUE)

pathToTheDESeq2Results = args[1]
signatureLength = args[2]
scoreMethod = args[3]
#scoreThreshold = args[4]
survivalType = args[4] 
outputFileName = args[5]
pathToHelperFile = args[6]
seed = args[7]
pathForBetaCoeff = args[8]

#Load the data
suppressPackageStartupMessages(library("DESeq2"))
load(pathToTheDESeq2Results)
expressionData = data.frame(assay(vsd))
#Extract the expression data for only protein coding genes
proteinCodingGenes = rownames(res.annot)[res.annot$biotype == "protein_coding"]
expressionData = expressionData[proteinCodingGenes, ]
print(dim(expressionData))
## Extract random gene signatrue of same length as our signature of interest
set.seed(seed)
select = sample(rownames(expressionData), signatureLength, replace = FALSE) 
## Extract expression data corresponding to this gene signature
expressionData = expressionData[select, ]
clinicalData = data.frame(colData(vsd))

#Make sure that the clinical data is in the same order as samples in the expression data
seIdx = match(rownames(clinicalData), colnames(expressionData))
expressionData = expressionData[, seIdx]

## Extract metadata 
### Modified from Dom's code in 3-survival v 0.4

suppressPackageStartupMessages(library("survival"))

# survival outcome 1: "d" for death and "ltrc" for left truncated and right censored
clinicalData$survival_d_ltrc = Surv(time  = as.numeric(difftime(clinicalData$DOE, clinicalData$DDiag))/365.25,
                                    time2 = as.numeric(difftime(clinicalData$DOC, clinicalData$DDiag))/365.25,
                                    event = ifelse(clinicalData$Dead == FALSE, 0, 1))

# survival outcome 2: "rd" for relapse or death, "ltrc" for left truncated and right censored
clinicalData$survival_rd_ltrc = Surv(time  = as.numeric(difftime(clinicalData$DOE, clinicalData$DDiag))/365.25,
                                     time2 = as.numeric(difftime(apply(clinicalData[, c("DOC", "DDistMets")],1,min,na.rm=TRUE), clinicalData$DDiag))/365.25,
                                     event = ifelse((clinicalData$Dead == FALSE & clinicalData$EventMet == "No"), 0, 1))

## Divide the samples into high-risk, low-risk groups based on gene expression and given (method, threshold) pair.
source(pathToHelperFile)
reducedSignature = applyMethod(expressionData, scoreMethod, pathToTheDESeq2Results, pathForBetaCoeff)
#expressionGroupsBasedOnSignature = as.factor(divideIntoCategories(reducedSignature, scoreThreshold))
clinicalData$Signature = reducedSignature
#clinicalData$SignatureGroups = expressionGroupsBasedOnSignature

## Perform survival analysis

#install.packages("survminer")
survival = ifelse((survivalType == "PFS"), "survival_rd_ltrc", "survival_d_ltrc")

## Multivariate analysis
covariates <- c("Sex + Age + as.character(Stage) + ECOG + treatment + as.numeric(Nclass)")
multiv_formulas <- sapply(covariates,
                          function(x) as.formula(paste(survival,'~ Signature +', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data.frame(clinicalData))})
# Extract data
multiv_results = list()
multiv_results <- lapply(multiv_models, saveSurvivalResults)
resMultiv <- t(as.data.frame(multiv_results, check.names = FALSE))
rownames(resMultiv) = c("Signature + Sex + Age + as.character(Stage) + ECOG + treatment + as.numeric(Nclass)")

## Save the results
write.table(resMultiv, file = outputFileName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
