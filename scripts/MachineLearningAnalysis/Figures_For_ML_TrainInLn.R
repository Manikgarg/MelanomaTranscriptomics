library("ggplot2")
library("dplyr")
library("openxlsx")

theme_set(ggpubr::theme_pubr(base_size=10, legend='bottom', x.text.angle = 45))

filePath = paste0("~/Desktop/Melanoma/TrainInLn/")
pattern = "resultTrain_"
tempFiles = list.files(path = filePath, pattern = pattern, recursive = TRUE)
currentResults = data.frame()
for (temp in tempFiles) {
  splitString = unlist(strsplit(tools::file_path_sans_ext(temp), split='_', fixed=TRUE))
  classifierName = splitString[length(splitString) -3]
  if (classifierName != "plr") {
    fullName = paste0(filePath, temp)
    tempResults = read.table(file = fullName, header = TRUE, quote = "", sep = "\t")
    tempResults = tempResults[, c("ROC" ,"ROCSD", "Sens", "SensSD", "Spec", "SpecSD")]
    tempResults$Resampling = splitString[length(splitString)]
    tempResults$PCA = splitString[length(splitString) -2]
    tempResults$ClassifierName = splitString[length(splitString) -4]
    tempResults$Signature = paste0(splitString[2:(length(splitString) -5)], collapse = "_")
    currentResults = rbind(currentResults, tempResults)
  }
}

currentResults = currentResults[currentResults$PCA!=TRUE, ]
currentResults$Signature_f = as.factor(currentResults$Signature)
levels(currentResults$Signature_f) = c("Clinical Covariates", "DecisionDx-\nMelanoma",
                                     "LMC_150", "Cam_121", "Cam_121+\nClinical Covariates")
currentResults$Signature_f <- factor(currentResults$Signature_f, 
                                     levels = c("Cam_121+\nClinical Covariates", 
                                                "Cam_121", "Clinical Covariates", 
                                                "LMC_150", "DecisionDx-\nMelanoma"))

currentResults$Resampling <- as.factor(currentResults$Resampling)
levels(currentResults$Resampling) <- c("Bootstrap", "10-Fold CV")
#write.xlsx(currentResults, "~/Desktop/Melanoma/githubUpload/Source_Data/Figs_3A_S7_S8_And_SupplementaryTable_S3A.xlsx", colNames = TRUE, rowNames = TRUE, append = TRUE)

#library("ggplot2")
linesDf<-currentResults%>%group_by(Signature_f, Resampling)%>%summarise(max=max(ROC), median=median(ROC), mean=mean(ROC))
#currentResults<-merge(currentResults, linesDf, by="Signature_f")
  
d1<- ggplot(currentResults, aes(x = ClassifierName, y = ROC, color = Resampling))+
  geom_point(size = 1, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = ROC - ROCSD, 
                    ymax = ROC + ROCSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual", name="Resampling\nmethod")+
  geom_hline(aes(yintercept=median, linetype=Signature_f, color=Resampling), data=linesDf)+
  scale_linetype(name="Median ROC")+
  ylim(c(0, 1))+
  ylab('AUROC')+
  xlab('Classifier name')+
  ggtitle("Training dataset = Lymph node")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)+
  ggpubr::theme_pubr(base_size=10, legend='right', x.text.angle = 90)
ggsave("~/Desktop/Melanoma/Figure3A'_ROC.png", device = "png", 
       width = 18.3, height = 10, units = "cm")

linesDf<-currentResults%>%group_by(Signature_f, Resampling)%>%summarise(max=max(Sens), median=median(Sens), mean=mean(Sens))

d2<- ggplot(currentResults, aes(x = ClassifierName, y = Sens, color = Resampling))+
  geom_point(size = 1, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = Sens - SensSD, 
                    ymax = Sens + SensSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  geom_hline(aes(yintercept=median, linetype=Signature_f, color=Resampling), data=linesDf)+
  scale_linetype(name="Median\nSensitivity")+
  ylim(c(0, 1))+
  ylab("Sensitivity")+
  xlab('Classifier name')+
  ggtitle("Training dataset = Lymph node")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)+
  ggpubr::theme_pubr(base_size=10, legend='right', x.text.angle = 90)
ggsave("~/Desktop/Melanoma/Figure3A'_Sens.png", device = "png", 
       width = 18.3, height = 10, units = "cm")

linesDf<-currentResults%>%group_by(Signature_f, Resampling)%>%summarise(max=max(Spec), median=median(Spec), mean=mean(Spec))
d3<- ggplot(currentResults, aes(x = ClassifierName, y = Spec, color=Resampling))+
  geom_point(size = 1, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = Spec - SpecSD, 
                    ymax = Spec + SpecSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  geom_hline(aes(yintercept=median, linetype=Signature_f, color=Resampling), data=linesDf)+
  scale_linetype(name="Median\nSpecificity")+
  ylim(c(0, 1))+
  ylab("Specificity")+
  xlab('Classifier name')+
  ggtitle("Training dataset = Lymph node")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)+
  ggpubr::theme_pubr(base_size=10, legend='right', x.text.angle = 90)
ggsave("~/Desktop/Melanoma/Figure3A'_Spec.png", device = "png", 
       width = 18.3, height = 10, units = "cm")

# Train p-values ----------------------------------------------------------
#### k-fold
tStatistics = data.frame()
for(metric in c("ROC", "Sens", "Spec")){
  tStatistic1 = t.test(x = currentResults[currentResults$Signature_f == "Cam_121+\nClinical Covariates", metric], y = currentResults[currentResults$Signature_f == "Clinical Covariates", metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic2 = t.test(x = currentResults[(currentResults$Signature_f == "Cam_121") & (currentResults$Resampling == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$Resampling == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic3 = t.test(x = currentResults[(currentResults$Signature_f == "DecisionDx-\nMelanoma") & (currentResults$Resampling == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$Resampling == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic4 = t.test(x = currentResults[(currentResults$Signature_f == "LMC_150") & (currentResults$Resampling == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$Resampling == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic5 = t.test(x = currentResults[(currentResults$Signature_f == "Cam_121") & (currentResults$Resampling == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "DecisionDx-\nMelanoma") & (currentResults$Resampling == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic6 = t.test(x = currentResults[(currentResults$Signature_f == "Cam_121") & (currentResults$Resampling == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "LMC_150") & (currentResults$Resampling == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")

  tStatistic = data.frame("Cam_121_Cov_Vs_Cov" = c(tStatistic1[["statistic"]][["t"]], tStatistic1[["p.value"]], tStatistic1[["conf.int"]], tStatistic1[["stderr"]]),
                          "Cam_121_Vs_Cov" = c(tStatistic2[["statistic"]][["t"]], tStatistic2[["p.value"]], tStatistic2[["conf.int"]], tStatistic2[["stderr"]]),
                          "DecisionDxMelanoma_Vs_Cov" = c(tStatistic3[["statistic"]][["t"]], tStatistic3[["p.value"]], tStatistic3[["conf.int"]], tStatistic3[["stderr"]]), 
                          "LMC_150_Vs_Cov" = c(tStatistic4[["statistic"]][["t"]], tStatistic4[["p.value"]], tStatistic4[["conf.int"]], tStatistic4[["stderr"]]),
                          "Cam_121_Vs_DecisionDxMelanoma" = c(tStatistic5[["statistic"]][["t"]], tStatistic5[["p.value"]], tStatistic5[["conf.int"]], tStatistic5[["stderr"]]),
                          "Cam_121_Vs_LMC_150" = c(tStatistic6[["statistic"]][["t"]], tStatistic6[["p.value"]], tStatistic6[["conf.int"]], tStatistic6[["stderr"]]))
  rownames(tStatistic) = paste(metric, c("t.statistic", "p.value", "conf.int1", "conf.int2", "SE"), sep = ".")
  tStatistics = rbind(tStatistics, tStatistic)
} 

write.table(tStatistics, "~/Desktop/Melanoma/train_kfold_PVal.tsv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

#### Bootstrap strap
tStatistics = data.frame()
for(metric in c("ROC", "Sens", "Spec")){
  tStatistic1 = t.test(x = currentResults[currentResults$Signature_f == "Cam_121+\nClinical Covariates", metric], y = currentResults[currentResults$Signature_f == "Clinical Covariates", metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic2 = t.test(x = currentResults[(currentResults$Signature_f == "Cam_121") & (currentResults$Resampling == "Bootstrap"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$Resampling == "Bootstrap"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic3 = t.test(x = currentResults[(currentResults$Signature_f == "DecisionDx-\nMelanoma") & (currentResults$Resampling == "Bootstrap"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$Resampling == "Bootstrap"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic4 = t.test(x = currentResults[(currentResults$Signature_f == "LMC_150") & (currentResults$Resampling == "Bootstrap"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$Resampling == "Bootstrap"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic5 = t.test(x = currentResults[(currentResults$Signature_f == "Cam_121") & (currentResults$Resampling == "Bootstrap"), metric], y = currentResults[(currentResults$Signature_f == "DecisionDx-\nMelanoma") & (currentResults$Resampling == "Bootstrap"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic6 = t.test(x = currentResults[(currentResults$Signature_f == "Cam_121") & (currentResults$Resampling == "Bootstrap"), metric], y = currentResults[(currentResults$Signature_f == "LMC_150") & (currentResults$Resampling == "Bootstrap"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  
  tStatistic = data.frame("Cam_121_Cov_Vs_Cov" = c(tStatistic1[["statistic"]][["t"]], tStatistic1[["p.value"]], tStatistic1[["conf.int"]], tStatistic1[["stderr"]]),
                          "Cam_121_Vs_Cov" = c(tStatistic2[["statistic"]][["t"]], tStatistic2[["p.value"]], tStatistic2[["conf.int"]], tStatistic2[["stderr"]]),
                          "DecisionDxMelanoma_Vs_Cov" = c(tStatistic3[["statistic"]][["t"]], tStatistic3[["p.value"]], tStatistic3[["conf.int"]], tStatistic3[["stderr"]]), 
                          "LMC_150_Vs_Cov" = c(tStatistic4[["statistic"]][["t"]], tStatistic4[["p.value"]], tStatistic4[["conf.int"]], tStatistic4[["stderr"]]),
                          "Cam_121_Vs_DecisionDxMelanoma" = c(tStatistic5[["statistic"]][["t"]], tStatistic5[["p.value"]], tStatistic5[["conf.int"]], tStatistic5[["stderr"]]),
                          "Cam_121_Vs_LMC_150" = c(tStatistic6[["statistic"]][["t"]], tStatistic6[["p.value"]], tStatistic6[["conf.int"]], tStatistic6[["stderr"]]))
  rownames(tStatistic) = paste(metric, c("t.statistic", "p.value", "conf.int1", "conf.int2", "SE"), sep = ".")
  tStatistics = rbind(tStatistics, tStatistic)
} 

write.table(tStatistics, "~/Desktop/Melanoma/train_Bootstrap_PVal.tsv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)