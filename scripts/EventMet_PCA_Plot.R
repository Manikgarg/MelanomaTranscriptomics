library("DESeq2")
library("ggplot2")
theme_set(ggpubr::theme_pubr(base_size=10, legend='bottom'))

load("./Downloads/de-7.Rdata")
expressionData = data.frame(assay(vsd))

clinicalData = data.frame(colData(vsd))
plotPCARiskGroup = function(geneExpressionDf, plotScreePlot, plotType){
  PCA = prcomp(geneExpressionDf)
  percentVar = round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  for (xComponent in 1:3){
    for (yComponent in (xComponent+1):4){
      dataGG = data.frame(xPComponent = PCA$x[,xComponent], yPComponent = PCA$x[,yComponent], 
                          Tissue_Code = clinicalData$Tissue_Code,
                          Metastases = clinicalData$EventMet)
      plot1 = ggplot(dataGG, aes(x = xPComponent, y = yPComponent, 
                                 shape = clinicalData$Tissue_Code, 
                                 color = clinicalData$EventMet)) +
        geom_point(alpha = 0.5)+
        labs(x = paste0("PC",xComponent," (", round(percentVar[xComponent],4),"%)"),
             y = paste0("PC",yComponent," (", round(percentVar[yComponent],4),"%)")) + 
        scale_color_manual(values = c("#1F78B4", "#E31A1C"), name = "Distant\nMetastases")+
        scale_shape(name = "Tissue")+
        coord_fixed()
      
      ggsave(filename = paste("./Desktop/Melanoma/", xComponent, "Vs", yComponent, "PlotForMetastases", 
                              plotType, ".png", sep = ""), width=8, units='cm')
      #dev.off()
    }}
  if(plotScreePlot & !file.exists(paste("./Desktop/Melanoma/ScreePlotForMetastases", plotType, ".pdf", sep = ""))){
    png(file=paste("./Desktop/Melanoma/ScreePlotForMetastases", plotType, ".png", sep = ""),
        width = 8, units = "cm")
    plot2 = screeplot(PCA, type = 'lines')
    dev.off()
    #Return both the plots as a list
    #Courtsey: https://stackoverflow.com/questions/35849703/returning-two-graphs-using-a-single-function-in-r
  }
}

plotPCARiskGroup(t(expressionData), plotScreePlot = TRUE, plotType = "AllGenes")

#Subset to top 1000 most variable genes

nTop = 1000
m = expressionData[, rownames(clinicalData[clinicalData$Tissue_Code=="Primary", ])]
Pvars <- rowVars(as.matrix(m))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(nTop, length(Pvars)))]
matSelected <- expressionData[select, ]
plotPCARiskGroup(t(matSelected), plotScreePlot = TRUE, plotType = "Top1000MostVariableGenes")

