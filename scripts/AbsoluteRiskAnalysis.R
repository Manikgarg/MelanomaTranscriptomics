data = read.table("./Downloads/LMCdata_Manik12feb2020.txt", sep = "\t", header = TRUE, quote = "")
rownames(data) = data$manik_id
data = data[, c("died_mm", "msstime", "ajccstage", "negativeScore121", "negativeScore70")]
data = data[!(rownames(data) == "4220077"), ]
data = data[complete.cases(data$ajccstage), ]

data$event = ifelse((data$died_mm == "Yes" & data$msstime < 5), "Yes", "No")
data = data[!(data$died_mm == "No" & data$msstime < 5), ]
table(data$event, data$ajccstage)

data$negativeScore121 = - data$negativeScore121
data$negativeScore70 = - data$negativeScore70

divideIntoCategories = function(reducedSignature, threshold){
  threshold = ifelse(threshold == "median", median(reducedSignature), as.numeric(quantile(reducedSignature,as.numeric(threshold))))
  categories = ifelse(reducedSignature >= threshold, "High", "Low")
  return(categories)
}

data$CutOff_median_121 = divideIntoCategories(data$negativeScore121, "median")
data$CutOff_0.33_121 = divideIntoCategories(data$negativeScore121, 0.33)
data$CutOff_0.25_121 = divideIntoCategories(data$negativeScore121, 0.25)
data$CutOff_0.25_0.75_121 = data$CutOff_0.25_121

data$CutOff_0.25_0.75_121[(data$negativeScore121 >= as.numeric(quantile(data$negativeScore121,as.numeric(0.25))))
                          & (data$negativeScore121 < as.numeric(quantile(data$negativeScore121,as.numeric(0.75))))] = "Mid"

write.table(data, "./Desktop/Melanoma/Data_For_Calculating_MSS_5y_AbsoluteRisk.tsv", sep = "\t", 
            col.names = TRUE, row.names = FALSE)

calculate5yearMSSAbsoluteRiskPerStage = function(data, threshold){
  absoluteRiskGroup = cbind(threshold, as.data.frame(table(GEP_class = data[[threshold]][data$event == "Yes"], Stage = data$ajccstage[data$event == "Yes"])))
  colnames(absoluteRiskGroup) = c("Threshold", "GEP_class", "Stage", "Event=Yes")
  totalGroup = as.data.frame(table(GEP_class = data[[threshold]], Stage = data$ajccstage))
  mergedGroups = cbind(absoluteRiskGroup, totalGroup)
  absoluteRiskGroup$Total = mergedGroups[, ncol(mergedGroups)]
  absoluteRiskGroup$Proportion = absoluteRiskGroup$`Event=Yes`/absoluteRiskGroup$Total
  return(absoluteRiskGroup)
}

absoluteRiskGroup = rbind(calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_median_121"), 
                          calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_0.33_121"), 
                          calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_0.25_121"), 
                          calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_0.25_0.75_121"))

write.table(absoluteRiskGroup, paste("./Desktop/Melanoma/MSS_5y_AbsoluteRisk.tsv", sep = ""), sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

require(survival)
require(survcomp)
survival_Object = Surv(data$msstime, as.numeric(factor(data$died_mm)))
pdf(file = "./Desktop/Melanoma/Leeds_MSS_KMcurve_0.33.pdf")
overall_fit = survfit(survival_Object ~ factor(data$CutOff_0.33_121))
plot(overall_fit, col = c("red", "blue"), xlab = "Time (years)", ylab = "Survival probability")
#pvalStr = res["Signature","p.value"]
#text(x= 5, y=0.1, paste0("p-value = ", pvalStr))
legend("topright", col = c("red", "blue"), legend = unique(data$CutOff_0.33_121), lty = 1)
dev.off()

sa<-summary(coxph(survival_Object ~ factor(data$CutOff_0.33_121)))
HR.OS<-c(P=sa$sctest[3],HR=sa$conf.int[1],CIlo=sa$conf.int[3],CIup=sa$conf.int[4])
print(HR.OS)

data$CutOff_0.33_70 = divideIntoCategories(data$negativeScore70, 0.33)
require(survival)
require(survcomp)
survival_Object = Surv(data$msstime, as.numeric(factor(data$died_mm)))
pdf(file = "./Desktop/Melanoma/Leeds_MSS_KMcurve_0.33_Cam70.pdf")
overall_fit = survfit(survival_Object ~ factor(data$CutOff_0.33_70))
plot(overall_fit, col = c("red", "blue"), xlab = "Time (years)", ylab = "Survival probability")
#pvalStr = res["Signature","p.value"]
#text(x= 5, y=0.1, paste0("p-value = ", pvalStr))
legend("topright", col = c("red", "blue"), legend = unique(data$CutOff_0.33_70), lty = 1)
dev.off()

sa<-summary(coxph(survival_Object ~ factor(data$CutOff_0.33_70)))
HR.OS<-c(P=sa$sctest[3],HR=sa$conf.int[1],CIlo=sa$conf.int[3],CIup=sa$conf.int[4])
print(HR.OS)

sa<-summary(coxph(survival_Object ~ data$negativeScore70))
HR.OS<-c(P=sa$sctest[3],HR=sa$conf.int[1],CIlo=sa$conf.int[3],CIup=sa$conf.int[4])
print(HR.OS)
