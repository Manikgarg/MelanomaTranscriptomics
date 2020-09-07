#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

pathToHelperFile = args[1]
pathTogeneToGroupMappings = args[2]
group = args[3]
pathToTheDESeq2Results = args[4]
pathToSurvivalMetadata = args[5]
scoreMethod = args[6]
scoreThreshold = args[7]
survivalType = args[8] 
outputFileName = args[9]
#outputPlotName = args[10]
pathForBetaCoeff = args[10]
pathForTilCounts = args[11]

source(pathToHelperFile)
geneToGroupMappings = read.table(pathTogeneToGroupMappings, sep = "\t", header = TRUE, quote = "")
data = extractExpressionAndClinicalDataSet(pathToTheDESeq2Results, geneToGroupMappings, group, pathToSurvivalMetadata)
expressionData = data[[1]]
clinicalData = data[[2]]

#Extract tilCounts
tilCount = read.table(pathForTilCounts, sep = "\t", header = TRUE, quote = "")
tilCount = tilCount[!duplicated(tilCount$R.Seq_sampleID), c("ClarkScore", "Scanned_file.me", "R.Seq_sampleID", "MIAScore")]
rownames(tilCount) = tilCount$R.Seq_sampleID
#Add tilCounts to Clinical data
rownames(clinicalData) = clinicalData$RNA.Seq.Sample
clinicalData$ClarkScore <- NA
clinicalData$ClarkScore[rownames(tilCount)%in%rownames(clinicalData)]= as.character(tilCount$ClarkScore[rownames(tilCount)%in%rownames(clinicalData)])
clinicalData$ClarkScore = as.factor(clinicalData$ClarkScore)

## Extract metadata 
### Modified from Dom's code in 3-survival v 0.4

suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("survminer"))

# survival outcome 1: "d" for death and "ltrc" for left truncated and right censored
clinicalData$survival_d_ltrc = Surv(time  = as.numeric(difftime(clinicalData$DOE, clinicalData$DDiag))/365.25,
                                    time2 = as.numeric(difftime(clinicalData$DOC, clinicalData$DDiag))/365.25,
                                    event = ifelse(clinicalData$Dead == FALSE, 0, 1))

# survival outcome 2: "rd" for relapse or death, "ltrc" for left truncated and right censored
clinicalData$survival_rd_ltrc = Surv(time  = as.numeric(difftime(clinicalData$DOE, clinicalData$DDiag))/365.25,
                                     time2 = as.numeric(difftime(apply(clinicalData[, c("DOC", "DDistMets")],1,min,na.rm=TRUE), clinicalData$DDiag))/365.25,
                                     event = ifelse((clinicalData$Dead == FALSE & clinicalData$EventMet == "No"), 0, 1))

## Divide the samples into high-risk, low-risk groups based on gene expression and given (method, threshold) pair.

reducedSignature = applyMethod(expressionData, scoreMethod, pathToTheDESeq2Results, pathForBetaCoeff)
expressionGroupsBasedOnSignature = as.factor(divideIntoCategories(reducedSignature, scoreThreshold))
clinicalData$Signature = reducedSignature
clinicalData$SignatureGroups = expressionGroupsBasedOnSignature

## Perform survival analysis

#install.packages("survminer")
survival = ifelse((survivalType == "PFS"), "survival_rd_ltrc", "survival_d_ltrc")

## Univariate analysis
## Inlcude univariate analysis on all of covariates as well
covariates <- c("Sex", "Age", "as.character(Stage)", "ECOG", "treatment", "as.numeric(Nclass)", "Signature", "ClarkScore")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste(survival,'~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = clinicalData)})
# Extract data
univ_results = list()
univ_results <- lapply(univ_models, saveSurvivalResults)
res <- t(as.data.frame(univ_results, check.names = FALSE))
rownames(res) = c("Sex", "Age", "as.character(Stage)", "ECOG", "Treatment", "as.numeric(Nclass)", "Signature", "TIL count")

# KM curve for signature
#pdf(file = outputPlotName)
#overall_fit = survfit(as.formula(paste(survival,'~ SignatureGroups')), data = clinicalData)
#plot(overall_fit, col = c("red", "blue"), xlab = "Time (years)", ylab = "Survival probability")
#pvalStr = survminer::surv_pvalue(fit = overall_fit, data = clinicalData)$pval
#text(x= 5, y=0.1, paste0("p-value = ", pvalStr))
#legend("topright", col = c("red", "blue"), legend = levels(expressionGroupsBasedOnSignature), lty = 1)
#dev.off()

## Multivariate analysis
covariates <- c("Sex", "Age", "as.character(Stage)", "ECOG", "treatment", "as.numeric(Nclass)", "Sex + Age + as.character(Stage) + ECOG + treatment + as.numeric(Nclass)", "Sex + Age + as.character(Stage) + ECOG + treatment + as.numeric(Nclass) + ClarkScore")
multiv_formulas <- sapply(covariates,
                          function(x) as.formula(paste(survival,'~ Signature +', x)))

multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data.frame(clinicalData))})
# Extract data
multiv_results = list()
multiv_results <- lapply(multiv_models, saveSurvivalResults)
resMultiv <- t(as.data.frame(multiv_results, check.names = FALSE))
#rownames(resMultiv) = c("Signature + Age", "Signature + as.character(Stage)", "Signature + Breslow", "Signature + ECOG", "Signature + Treatment", "Signature + as.numeric(Nclass)", "Signature + Ulc", "Signature + as.character(Stage) + Breslow + ECOG + Treatment", "Signature + Age + as.character(Stage) + Breslow + ECOG + Treatment + as.numeric(Nclass)", "Signature + Age + as.character(Stage) + Breslow + ECOG + Treatment + as.numeric(Nclass) + Ulc", "Signature + Age + Breslow + ECOG + Treatment + as.numeric(Nclass) + Ulc")
rownames(resMultiv) = c("Signature + Sex", "Signature + Age", "Signature + as.character(Stage)", "Signature + ECOG", "Signature + treatment", "Signature + as.numeric(Nclass)", "Signature + Sex + Age + as.character(Stage) + ECOG + treatment + as.numeric(Nclass)", "Signature + Sex + Age + as.character(Stage) + ECOG + treatment + as.numeric(Nclass) + TIL count")
res <- rbind(res, resMultiv)

## Save the results
write.table(res, file = outputFileName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
