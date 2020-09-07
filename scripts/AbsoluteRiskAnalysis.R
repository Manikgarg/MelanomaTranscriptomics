library("openxlsx")

data = read.table("~/Downloads/LMCdata_Manik12feb2020.txt", sep = "\t", header = TRUE, quote = "")
rownames(data) = data$manik_id
data = data[, c("died_mm", "msstime", "ajccstage", "negativeScore121")] #, "negativeScore70")]
data = data[!(rownames(data) == "4220077"), ]
data = data[complete.cases(data$ajccstage), ]

data$event = ifelse((data$died_mm == "Yes" & data$msstime < 5), "Yes", "No")
data = data[!(data$died_mm == "No" & data$msstime < 5), ]
table(data$event, data$ajccstage)

data$negativeScore121 = - data$negativeScore121
#data$negativeScore70 = - data$negativeScore70

divideIntoCategories = function(reducedSignature, threshold){
  threshold = ifelse(threshold == "median", median(reducedSignature), as.numeric(quantile(reducedSignature,as.numeric(threshold))))
  categories = ifelse(reducedSignature >= threshold, "High", "Low")
  return(categories)
}

#data$CutOff_median_121 = divideIntoCategories(data$negativeScore121, "median")
data$CutOff_0.33_121 = divideIntoCategories(data$negativeScore121, 0.33)
#data$CutOff_0.25_121 = divideIntoCategories(data$negativeScore121, 0.25)
#data$CutOff_0.25_0.75_121 = data$CutOff_0.25_121

#data$CutOff_0.25_0.75_121[(data$negativeScore121 >= as.numeric(quantile(data$negativeScore121,as.numeric(0.25))))
#                          & (data$negativeScore121 < as.numeric(quantile(data$negativeScore121,as.numeric(0.75))))] = "Mid"

write.table(data, "./Desktop/Melanoma/Data_For_Calculating_MSS_5y_AbsoluteRisk.tsv", sep = "\t", 
            col.names = TRUE, row.names = FALSE)

write.xlsx(data, "~/Desktop/Melanoma/githubUpload/Source_Data/Table_1.xlsx", colNames = TRUE, rowNames = TRUE, append = TRUE)

calculate5yearMSSAbsoluteRiskPerStage = function(data, threshold){
  absoluteRiskGroup = cbind(threshold, as.data.frame(table(GEP_class = data[[threshold]][data$event == "Yes"], Stage = data$ajccstage[data$event == "Yes"])))
  colnames(absoluteRiskGroup) = c("Threshold", "GEP_class", "Stage", "Event=Yes")
  totalGroup = as.data.frame(table(GEP_class = data[[threshold]], Stage = data$ajccstage))
  mergedGroups = cbind(absoluteRiskGroup, totalGroup)
  absoluteRiskGroup$Total = mergedGroups[, ncol(mergedGroups)]
  absoluteRiskGroup$Proportion = absoluteRiskGroup$`Event=Yes`/absoluteRiskGroup$Total
  return(absoluteRiskGroup)
}

absoluteRiskGroup = rbind(#calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_median_121"), 
                          calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_0.33_121")#, 
                          #calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_0.25_121"), 
                          #calculate5yearMSSAbsoluteRiskPerStage(data, "CutOff_0.25_0.75_121")
                          )

write.table(absoluteRiskGroup, paste("~/Desktop/Melanoma/MSS_5y_AbsoluteRisk.tsv", sep = ""), sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
