library("ggplot2")
library("openxlsx")

theme_set(ggpubr::theme_pubr(base_size=10, legend='bottom', x.text.angle = 45))

filePath = paste0("~/Desktop/Melanoma/kFold_1000repeats/")
pattern = "resultTrain_"
tempFiles = list.files(path = filePath, pattern = pattern, recursive = TRUE)
currentResults = data.frame()
for (temp in tempFiles) {
  splitString = unlist(strsplit(temp, split='_', fixed=TRUE))
  classifierName = splitString[length(splitString) -3]
  if (classifierName != "plr") {
    fullName = paste0(filePath, temp)
    tempResults = read.table(file = fullName, header = TRUE, quote = "", sep = "\t")
    tempResults = tempResults[, c("ROC" ,"ROCSD", "Sens", "SensSD", "Spec", "SpecSD")]
    tempResults$PCA = splitString[length(splitString) -1]
    tempResults$ClassifierName = splitString[length(splitString) -3]
    tempResults$Signature = paste0(splitString[2:(length(splitString) -4)], collapse = "_")
    tempResults$CV = "10-Fold CV"
    currentResults = rbind(currentResults, tempResults)
  }
}

filePath = paste0("~/Desktop/Melanoma/loocv/")
pattern = "resultTrain_"
tempFiles = list.files(path = filePath, pattern = pattern, recursive = TRUE)
#currentResults = data.frame()
for (temp in tempFiles) {
  splitString = unlist(strsplit(temp, split='_', fixed=TRUE))
  classifierName = splitString[length(splitString) -3]
  if (classifierName != "plr") {
    fullName = paste0(filePath, temp)
    tempResults = read.table(file = fullName, header = TRUE, quote = "", sep = "\t")
    tempResults = tempResults[, c("ROC" ,"Sens", "Spec")]
    tempResults$ROCSD = 0
    tempResults$SensSD = 0
    tempResults$SpecSD = 0
    tempResults$PCA = splitString[length(splitString) -1]
    tempResults$ClassifierName = splitString[length(splitString) -3]
    tempResults$Signature = paste0(splitString[3:(length(splitString) -4)], collapse = "_")
    tempResults$CV = "LOOCV"
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

write.xlsx(currentResults, "~/Desktop/Melanoma/githubUpload/Source_Data/Figs_3A_S7_S8_And_SupplementaryTable_S3A.xlsx", colNames = TRUE, rowNames = TRUE, append = TRUE)

#library("ggplot2")
d1<- ggplot(currentResults[(currentResults$Signature_f%in%c("Cam_121+\nClinical Covariates", 
                                               "Cam_121", "Clinical Covariates")), ], aes(x = ClassifierName, y = ROC, color = CV))+
  geom_point(size = 2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = ROC - ROCSD, 
                    ymax = ROC + ROCSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  ylim(c(0, 1))+
  xlab("Classifier name")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)
  #ggpubr::theme_pubr(base_size=10, legend='bottom')
ggsave("~/Desktop/Melanoma/kFold_Vs_LOOCV_mainOnly.png", device = "png", 
       width = 16, height = 10, units = "cm")

d0<- ggplot(currentResults[currentResults$CV == "10-Fold CV", ], 
            aes(x = ClassifierName, y = ROC, color = ClassifierName))+
  geom_point(size = 2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = ROC - ROCSD, 
                    ymax = ROC + ROCSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  ylim(c(0, 1))+
  ylab("AUROC")+
  xlab("Classifier name")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)
ggsave("~/Desktop/Melanoma/kFold.png", device = "png", 
       width = 16, height = 10, units = "cm")

d0_1<- ggplot(currentResults[(currentResults$CV == "10-Fold CV")&
                               (currentResults$Signature_f%in%c("Cam_121+\nClinical Covariates", 
                                                                "Cam_121", "Clinical Covariates")), ], 
            aes(x = ClassifierName, y = ROC, color = ClassifierName))+
  geom_point(size = 2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = ROC - ROCSD, 
                    ymax = ROC + ROCSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual", name="Classifier name")+
  ylim(c(0, 1))+
  ylab("AUROC")+
  xlab("Classifier name")+
  ggtitle("Training Dataset = Primary melanoma")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)
ggsave("~/Desktop/Melanoma/Figure3A.png", device = "png", 
       width = 12, height = 10, units = "cm")

#library("ggplot2")
d2<- ggplot(currentResults, aes(x = ClassifierName, y = Sens, color = CV))+
  geom_point(size = 2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = Sens - SensSD, 
                    ymax = Sens + SensSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  ylim(c(0, 1))+
  xlab("Classifier name")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)
  #ggpubr::theme_pubr(base_size=10, legend='bottom')

d02 <- ggplot(currentResults[(currentResults$CV == "10-Fold CV") & 
                               (currentResults$Signature_f%in%c("Cam_121+\nClinical Covariates",
                                                                "Cam_121", "Clinical Covariates")), ], 
              aes(x = ClassifierName, y = Sens, color = ClassifierName))+
  geom_point(size = 2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = Sens - SensSD, 
                    ymax = Sens + SensSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  ylim(c(0, 1))+
  ylab("Sensitivity")+
  xlab("Classifier name")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)
ggsave("~/Desktop/Melanoma/kFoldSens.png", device = "png", 
       width = 16, height = 10, units = "cm")

#library("ggplot2")
d3<- ggplot(currentResults, aes(x = ClassifierName, y = Spec, color = CV))+
  geom_point(size = 2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = Spec - SpecSD, 
                    ymax = Spec + SpecSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  ylim(c(0, 1))+
  xlab("Classifier name")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)
  #ggpubr::theme_pubr(base_size=10, legend='bottom')

d03 <- ggplot(currentResults[(currentResults$CV == "10-Fold CV") & 
                               (currentResults$Signature_f%in%c("Cam_121+\nClinical Covariates",
                                                                "Cam_121", "Clinical Covariates")), ], 
              aes(x = ClassifierName, y = Spec, color = ClassifierName))+
  geom_point(size = 2, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin = Spec - SpecSD, 
                    ymax = Spec + SpecSD), 
                width = .1,
                position=position_dodge(width=0.5)) +
  scale_color_brewer(palette = "Set2", type = "qual")+
  ylim(c(0, 1))+
  ylab("Specificity")+
  xlab("Classifier name")+
  #theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_grid(.~Signature_f, scales = "free_x", drop = TRUE)
ggsave("~/Desktop/Melanoma/kFoldSpec.png", device = "png", 
       width = 16, height = 10, units = "cm")

mainResults = data.frame()
## Select the hyperparameters giving maximum value for each signature
for (i in unique(currentResults$Signature)) {
  currentResultsFori = currentResults[currentResults$Signature == i, ]
  idx = which.max(currentResultsFori[, "ROC"])
  mainResults = rbind(mainResults, currentResultsFori[idx, ])
}

# PredictionOnTest --------------------------------------------------------

extractTestPred <- function(group, filePath, mainResults){
  currentResultsFori = mainResults[mainResults$Signature == group, ]
  if(mainResults$CV[mainResults$Signature == group] == "LOOCV"){
    fileName = paste("predTest_loocv", currentResultsFori$Signature, 
                     currentResultsFori$ClassifierName, "pca", 
                     currentResultsFori$PCA, "3.tsv", sep = "_")
  }else{
    fileName = paste("predTest", currentResultsFori$Signature, 
                     currentResultsFori$ClassifierName, "pca", 
                     currentResultsFori$PCA, "3.tsv", sep = "_")
  }
  predTest = read.table(paste0(filePath, fileName), sep = "\t", header = TRUE, quote = "")
  predTest$Signature = mainResults$Signature_f[mainResults$Signature == group]
  return(predTest)
}

library("pROC")
predTest = data.frame()
for (i in unique(mainResults$Signature)) {
  if(mainResults$CV[mainResults$Signature == i] == "LOOCV"){
    filePath = paste0("~/Desktop/Melanoma/loocv/")
  }else{
    filePath = paste0("~/Desktop/Melanoma/kFold_1000repeats/")
  }
  temp = extractTestPred(i, filePath, mainResults)
  temp$roc = round(auc(roc(temp$obs, temp$Yes)), 4)
  temp$Signature_ROC = paste(temp$Signature, "\n(AUROC=", round(temp$roc, 2), ")", sep = "")
  temp$CV = mainResults$CV[mainResults$Signature == i]
  temp$ClassifierName = mainResults$ClassifierName[mainResults$Signature == i]
  predTest = rbind(predTest, temp)
}

#predTest$Signature_ROC = as.factor(predTest$Signature_ROC)
predTest$Signature_ROC <- factor(predTest$Signature_ROC, 
                                     levels = unique(predTest$Signature_ROC[order(predTest$roc, decreasing = TRUE)]))

write.xlsx(predTest, "~/Desktop/Melanoma/githubUpload/Source_Data/Figs_3B_3C_3D_And_SupplementaryTable_S3B.xlsx", colNames = TRUE, rowNames = TRUE, append = TRUE)

#library("ggplot2")
library("plotROC")

g <- ggplot(predTest[predTest$Signature%in%c("Cam_121+\nClinical Covariates", "Cam_121", "Clinical Covariates"), ],          
  aes(m=Yes, d=factor(obs, levels = c("No", "Yes")), color = Signature_ROC)) +
  geom_roc(hjust = -0.4, vjust = 1.5, n.cuts=0) +
  #coord_equal() +
  style_roc(xlab = "False Positive Rate (1-Specificity)", ylab = "True Positive Rate (Sensitivity)") +
  geom_abline(intercept = 0, slope = 1, col = "gray60", lty = 2)+
  scale_color_brewer(type = "qual", palette = "Paired", name = "Signature") +
  #scale_linetype_manual(values = c(3, 1), name = "Classifier")+
  ggpubr::theme_pubr(base_size=10, legend = "right", x.text.angle = 45) +
  theme(legend.position=c(.825,.325), legend.text = element_text(size=6), 
        legend.title = element_text(size = 6))+
  ggtitle("Validation Dataset = Lymph node")+
  coord_fixed()

g

ggsave("~/Desktop/Melanoma/MLLnValidation.png", device = "png", units='cm', width = 8)

ggpubr::ggarrange(d1, g, ncol = 1, nrow = 2, 
                  labels = c("A", "B"),
                  width = c(2, 1))
ggsave("~/Desktop/Melanoma/MLkFold.pdf", device = "pdf", units='cm', width = 16)

ggpubr::ggarrange(d2, d3, ncol = 1, nrow = 2, 
                  labels = c("A", "B"), align = "v", common.legend = TRUE)
ggsave("~/Desktop/Melanoma/MLkFoldSensSpec.pdf", device = "pdf", units='cm', 
       width = 16)

#g + annotate("text", x=c(0.75, 0.75, 0.75, 0.75), y=c(0.20, 0.25, 0.30, 0.35), label=paste("AUC =", round((calc_auc(g))$AUC, 4)))


# See if the ROC curves are stistically different from each oth --------

pVals = data.frame()
for (i in unique(predTest$Signature[predTest$Signature!="Clinical Covariates"])) {
  for (j in unique(predTest$Signature[predTest$Signature!=i])) {
    rocResults = roc.test(roc1 = roc(predTest$obs[predTest$Signature==i], 
                                     predTest$Yes[predTest$Signature==i]), 
                          roc2 = roc(predTest$obs[predTest$Signature==j], 
                                     predTest$Yes[predTest$Signature==j]), 
                          alternative = "greater")
    
    temp = data.frame("Comparison"= paste(i, ">", j, sep = " "), 
                      "p-val" = rocResults[["p.value"]],
                      "statistic" = rocResults[["statistic"]])
    pVals = rbind(temp, pVals) 
  }
}

write.table(pVals, "~/Desktop/Melanoma/pValues.tsv", sep = "\t", col.names = TRUE, 
            row.names = FALSE, quote = FALSE)


# Comparison between Cam_121 and Clinical Covariates test labels ----------

predTest$predLabel = ifelse(predTest$Yes>0.5, "Yes", "No")

combinedPredictions = data.frame("Cam_121" = predTest$predLabel[predTest$Signature=="Cam_121"], 
                                 "CCov" = predTest$predLabel[predTest$Signature == "Clinical Covariates"], 
                                 "TrueLabel" = predTest$obs[predTest$Signature=="Cam_121"])

##Compare the observations for each sample
for (i in 1:nrow(combinedPredictions)) {
  if (combinedPredictions$Cam_121[i] == combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] != combinedPredictions$TrueLabel[i])
  {
    combinedPredictions$Comparison[i] = "E"
  }
  else if (combinedPredictions$Cam_121[i] != combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] == combinedPredictions$TrueLabel[i]) 
  {
    combinedPredictions$Comparison[i] = "C"
  }
  else if (combinedPredictions$Cam_121[i] != combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] != combinedPredictions$TrueLabel[i]) 
  {
    combinedPredictions$Comparison[i] = "None"
  }
  else if(combinedPredictions$Cam_121[i] == combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] == combinedPredictions$TrueLabel[i]) 
  {
    combinedPredictions$Comparison[i] = "CE"
  }
}

library("VennDiagram")
x = list(A = rownames(combinedPredictions)[combinedPredictions$Comparison %in% c("C", "CE")], 
         B = rownames(combinedPredictions)[combinedPredictions$Comparison %in% c("E", "CE")])

png("Cam_121_Vs_Cov.png", width = 8, units = "cm")
p = venn.diagram(x, cex = 2, cat.cex = 0, fill = c("mistyrose", "lightskyblue1"), col = c("mistyrose", "lightskyblue1"), alpha = c(0.5, 0.5), filename = NULL)
grid.draw(p)
dev.off()

# Comparison between Cam_121 + Clinical Covariates and Clinical Covariates test labels ----------

combinedPredictions = data.frame("Cam_121" = predTest$predLabel[predTest$Signature=="Cam_121+\nClinical Covariates"], 
                                 "CCov" = predTest$predLabel[predTest$Signature == "Clinical Covariates"], 
                                 "TrueLabel" = predTest$obs[predTest$Signature=="Cam_121+\nClinical Covariates"])

##Compare the observations for each sample
for (i in 1:nrow(combinedPredictions)) {
  if (combinedPredictions$Cam_121[i] == combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] != combinedPredictions$TrueLabel[i])
  {
    combinedPredictions$Comparison[i] = "E"
  }
  else if (combinedPredictions$Cam_121[i] != combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] == combinedPredictions$TrueLabel[i]) 
  {
    combinedPredictions$Comparison[i] = "C"
  }
  else if (combinedPredictions$Cam_121[i] != combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] != combinedPredictions$TrueLabel[i]) 
  {
    combinedPredictions$Comparison[i] = "None"
  }
  else if(combinedPredictions$Cam_121[i] == combinedPredictions$TrueLabel[i] & combinedPredictions$CCov[i] == combinedPredictions$TrueLabel[i]) 
  {
    combinedPredictions$Comparison[i] = "CE"
  }
}

library("VennDiagram")
x = list(A = rownames(combinedPredictions)[combinedPredictions$Comparison %in% c("C", "CE")], 
         B = rownames(combinedPredictions)[combinedPredictions$Comparison %in% c("E", "CE")])

png("Cam_121_Cov_Vs_Cov.png", width = 8, units = "cm")
p = venn.diagram(x, cex = 2, cat.cex = 0, fill = c("mistyrose", "lightskyblue1"), col = c("mistyrose", "lightskyblue1"), alpha = c(0.5, 0.5), filename = NULL)
grid.draw(p)
dev.off()


# Train p-values ----------------------------------------------------------

tStatistics = data.frame()
for(metric in c("ROC", "Sens", "Spec")){
  tStatistic1 = t.test(x = currentResults[currentResults$Signature_f == "Cam_121+\nClinical Covariates", metric], y = currentResults[currentResults$Signature_f == "Clinical Covariates", metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic2 = t.test(x = currentResults[(currentResults$Signature_f == "Cam_121") & (currentResults$CV == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$CV == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic3 = t.test(x = currentResults[(currentResults$Signature_f == "DecisionDx-\nMelanoma") & (currentResults$CV == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$CV == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  tStatistic4 = t.test(x = currentResults[(currentResults$Signature_f == "LMC_150") & (currentResults$CV == "10-Fold CV"), metric], y = currentResults[(currentResults$Signature_f == "Clinical Covariates") & (currentResults$CV == "10-Fold CV"), metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  
  #tStatistic4 = t.test(x = currentResults[currentResults$Signature == "Signature_overlap_DASLarray_genes", metric], y = currentResults[currentResults$Signature == "Gerami_genes" & currentResults$ClassifierName == "rf", metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  #tStatistic5 = t.test(x = currentResults[currentResults$Signature == "Signature_overlap_DASLarray_genes", metric], y = currentResults[currentResults$Signature == "LMC_150_genes" & currentResults$ClassifierName == "glmnet", metric], var.equal = FALSE, paired = FALSE, alternative = "greater")
  
  tStatistic = data.frame("Cam_121_Cov_Vs_Cov" = c(tStatistic1[["statistic"]][["t"]], tStatistic1[["p.value"]], tStatistic1[["conf.int"]], tStatistic1[["stderr"]]),
                          "Cam_121_Vs_Cov" = c(tStatistic2[["statistic"]][["t"]], tStatistic2[["p.value"]], tStatistic2[["conf.int"]], tStatistic2[["stderr"]]),
                          "DecisionDxMelanoma_Vs_Cov" = c(tStatistic3[["statistic"]][["t"]], tStatistic3[["p.value"]], tStatistic3[["conf.int"]], tStatistic3[["stderr"]]), 
                          "LMC_150_Vs_Cov" = c(tStatistic4[["statistic"]][["t"]], tStatistic4[["p.value"]], tStatistic4[["conf.int"]], tStatistic4[["stderr"]])
                          #"Cam_121_Cov_Vs_DecisionDxMelanoma" = c(tStatistic4[["p.value"]], tStatistic4[["conf.int"]], tStatistic4[["SEerr"]]),
                          #"Cam_121_Cov_Vs_LMC_150" = c(tStatistic5[["p.value"]], tStatistic5[["conf.int"]], tStatistic5[["SEerr"]])
  )
  rownames(tStatistic) = paste(metric, c("t.statistic", "p.value", "conf.int1", "conf.int2", "SE"), sep = ".")
  tStatistics = rbind(tStatistics, tStatistic)
} 

write.table(tStatistics, "~/Desktop/Melanoma/trainKFoldPVal.tsv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)


# Feature Importance Score ------------------------------------------------

library(ggplot2)
load("~/Desktop/Melanoma/kFold_1000repeats/trainedModel_Signature_overlap_DASLarray_genes_FALSE_rf_.Rdata")
fic = data.frame(resultsWithTrainedModels[["model"]][["finalModel"]][["importance"]])

write.xlsx(fic, "~/Desktop/Melanoma/githubUpload/Source_Data/Fig_S10.xlsx", colNames = TRUE, rowNames = TRUE, append = TRUE)
# 
fic$Gene = as.factor(rownames(fic))

fic$Gene.f = factor(fic$Gene, levels = levels(fic$Gene))

p <- ggplot(fic, aes(x= Gene, y=MeanDecreaseGini)) + 
  geom_point() + 
  scale_x_discrete(limits=c("KRTAP19.6","NXPH1","HNF4A","OR5AK2","NPPC","HTR3B","AIM2","KLHDC8A","CCR3",
                            "UTS2","SLC17A3","CTCFL","ZNF560","MTUS2","JPH3","NTRK1","CAMK2B","LRRC31",
                            "ALK","FREM2","ANGPTL7","PRF1","ABCC12","TDRD12","RAB15","LGI4","HHATL","HEPACAM",
                            "CACNG4","RP1","KNDC1","WIF1","ART1","DIO3OS","LAG3","SIGLEC12","FLOT1","COL24A1",
                            "PLA2G2D","DPEP3","GABBR2","GLDC","CTAGE1","GPR39","SFMBT2","KCNT1","ZIM2","DYSF",
                            "TFF2","C8G","GPR63","GFRA1","CPN2","MCF2","RGS7","ANK1","EEF1A2","SLCO5A1","CETP",
                            "KIF19","CASQ1","TTN","BMX","SLIT1","SOCS1","CCR5","DOCK3","F7","WDR49","KYNU",
                            "KSR2","GBP5","FASLG","PLCZ1","ERGIC2","RPL7L1","MFSD6L","ACOT1","AGXT","ME3",
                            "PROCA1","DMAP1","PHEX","TLR4","VAMP5","MIR222","TRNAU1AP","GRAMD1B","MLF1",
                            "JARID2","GCH1","ZNF697","OR52K1","CDC5L","KLC4","EHBP1L1","TAS2R60","SPAG6",
                            "TNFRSF10A","CLIC5","BIN1","PTPRG","AGBL4","RNF213","FNBP1L","TRIM22","NFYA",
                            "CENPQ","TMCC3","SLC29A1","ARRDC1","SOX4","RBBP4","PIK3R6","FLVCR2","TUBB",
                            "P2RY14","OCIAD2","CUL7","PARP11","STMN1", "ECOG.Symptomatic","Stage.IIIC","Stage.IIIA","Bres..2.4mm","Bres....2.0.mm",
                            "Stage.IIIB","Stage.IIC","treatment.Bevacizumab", "Stage.IIB",
                            "Bres..4.0mm"),
                   labels=c("Bres..4.0mm" = "Breslow > 4mm", 
                            "Bres..2.4mm" = "Breslow > 2mm-4mm",
                            "Bres....2.0.mm"= "Breslow <= 2mm",
                            "treatment.Bevacizumab" = "Treatment = Bevacizumab",
                            "ECOG.Symptomatic" = "ECOG = Symptomatic",
                            "Stage.IIIC" = "Stage IIIC",
                            "Stage.IIIA" = "Stage IIIA",
                            "Stage.IIIB" = "Stage IIIB",
                            "Stage.IIC" = "Stage IIC"
                            )) +
  xlab("Cam_121 + Clinical Covariates") + 
  ylab("Feature Importance Score (Mean Decrease Gini)") +
  coord_flip()+
  #theme_bw()+
  theme(axis.text.y = element_text(size = 7, face = "italic", angle = 0),
        axis.text.x = element_text(size = 7, angle = 0),
        legend.position="none") 
#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
#stat_summary(fun.data="mean_sdl", 
#               geom="crossbar", width=0.2 )
ggsave("~/Desktop/Melanoma/Cam_121_Cov_FeatureImportanceScore.png", device = "png", width=16)

