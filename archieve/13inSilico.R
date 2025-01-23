#library(nnls)
library(pheatmap)
library(RColorBrewer)
library(CIBERSORT)
#abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/all.data.RData")
load(file = "./dataset/COVID/RF_meanGini.rds")
load(file = "./dataset/COVID/RF_meanGini2.rds")
load(file = "./dataset/COVID/SVM_sigfeature.rds")
niteration <- 1220

load(file = "./dataset/COVID/ABIS_RF_sigMat_adjPMN.rds")
load(file = "./dataset/COVID/all.adjusted_PMN.RData")
all.tpm <- adjusted.tpm
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
load(file = "./dataset/COVID/RF_meanGini_adjPMN.rds")
load(file = "./dataset/COVID/RF_meanGini2_adjPMN.rds")
load(file = "./dataset/COVID/topVarGenes.rds")



load(file = "./dataset/COVID/all_adjusted.RData")
load(file = "./dataset/COVID/RF_meanGini_adj.rds")
load( file = "./dataset/COVID/SigMatrix_adj.rds")
all.tpm <- adjusted.tpm
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
niteration <- 2280
load(file = "./dataset/COVID/pal.celltype.rds")

### RF meanGini signatures####
for (i in c(100, 250, 500, 1000, 2000, 3000)) {
  keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:i],]), rownames(all.tpm))
  tpm.decon <- all.tpm[keep,]
  cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:i],], as.matrix(tpm.decon), QN = F, perm = 100)
  cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))
  
  pal.celltype <- brewer.pal(8, "Set1")
  pal.celltype <- colorRampPalette(pal.celltype)(29)  
  colnames(cibersort.df) <- all.pdata$id
  mat_col <- data.frame(group = all.pdata$cell)
  rownames(mat_col) <- colnames(cibersort.df)
  mat_colors <- list(group = pal.celltype)
  names(mat_colors$group) <- unique(all.pdata$cell)
  #cibersort.df <- round(cibersort.df*100)
  pheatmap(cibersort.df[order(rownames(cibersort.df)),order(all.pdata$cell)], 
           cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
           show_colnames     = F,show_rownames     = T, border_color = "black", 
           annotation_col    = mat_col,annotation_colors = mat_colors,
           main = paste("RF meanGini", i , sep = " ") )
  print(paste(i, "Iteration Done", sep = " "))
}

#from iteration results, #350 is the best performance
niteration <- 2280
keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))


colnames(cibersort.df) <- all.pdata$id
mat_col <- data.frame(group = all.pdata$cell)
rownames(mat_col) <- colnames(cibersort.df)
mat_colors <- list(group = pal.celltype)
#names(mat_colors$group) <- unique(all.pdata$cell)
#cibersort.df <- round(cibersort.df*100)
pheatmap(cibersort.df[,order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = paste("RF meanGini", niteration , sep = " ") )



niteration <- 500
keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))
pheatmap(cibersort.df[,order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = paste("RF meanGini", niteration , sep = " ") )

niteration <- 1000
keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon), QN= F)
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))
pheatmap(cibersort.df[,order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = paste("RF meanGini", niteration , sep = " ") )

niteration <- 2000
keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon), QN= F)
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))
pheatmap(cibersort.df[,order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = paste("RF meanGini", niteration , sep = " ") )

niteration <- 3000
keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon), QN= F)
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))
pheatmap(cibersort.df[,order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = paste("RF meanGini", niteration , sep = " ") )






##benchmarking ABIS #####
dim(all.tpm)
write.table(all.tpm, file = "./dataset/COVID/all.tpm_refMat.txt", sep = "\t")
#install.packages(c("shiny", "MASS", "preprocessCore"), dependencies = TRUE)
shiny::runGitHub("ABIS", user="giannimonaco")

decon.abis <- read.delim("./dataset/COVID/all.tpm_refMat.txt_deconvolution_ABIS.txt", sep = "\t")
decon.abis[, 1:10]
decon.abis[decon.abis<0]  <- 0
decon.abis <- sweep(decon.abis,2,colSums(decon.abis),`/`)
decon.abis <- round(decon.abis, 4)
decon.abis <- decon.abis[, grep("GSM", all.pdata$id)]
rownames(decon.abis) <- c("Mono_C", "NK", "CD8Tm", "CD4Tnaive", "CD8Tnaive", "Bnaive", "CD4Tm", "MAIT", "Tgd1", "PMN", "Tgd2", "Basophils", "Mono_NC", "Bm", "mDC", "pDC", "Plasmablasts")
decon.abis <- rbind(decon.abis, colSums(decon.abis[c("Tgd1", "Tgd2"),]) )
rownames(decon.abis)[18] <- "Tgd"
decon.abis <- decon.abis[c(1:8,10,12:18),]

decon.df <- cibersort.df[,grep("GSM", colnames(cibersort.df))]

pdata.benchmark <- all.pdata[grep("GSM", colnames(cibersort.df)),]
pdata.benchmark <- pdata.benchmark[grep("CD8Teff|Th|Tfh|Treg|CD4Teff|Bex", pdata.benchmark$cell, invert = T),]
pdata.benchmark$cell <- gsub("CD8Tcm", "CD8Tm", pdata.benchmark$cell)
pdata.benchmark$cell <- gsub("CD8Tem", "CD8Tm", pdata.benchmark$cell)
pdata.benchmark$cell <- gsub("VD2", "Tgd", pdata.benchmark$cell)
pdata.benchmark$cell <- gsub("Bnsm", "Bm", pdata.benchmark$cell)
pdata.benchmark$cell <- gsub("Bsm", "Bm", pdata.benchmark$cell)
pdata.benchmark$cell <- gsub("Mono_I", "Mono_NC", pdata.benchmark$cell)

keep <- intersect(rownames(pdata.benchmark), colnames(decon.abis))
decon.abis <- decon.abis[,keep]
keep <- intersect(rownames(pdata.benchmark), colnames(decon.df))
decon.df <- decon.df[,keep]
unique(pdata.benchmark$cell)
unique(rownames(decon.abis))
decon.abis <- decon.abis[-7,]

pal.celltype <- brewer.pal(8, "Accent")
pal.celltype <- colorRampPalette(pal.celltype)(15)  
colnames(decon.abis) <- pdata.benchmark$id
mat_col <- data.frame(group = pdata.benchmark$cell)
rownames(mat_col) <- colnames(decon.abis)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata.benchmark$cell)
pheatmap(decon.abis[order(rownames(decon.abis)),order(pdata.benchmark$cell) ],
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = NA, 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "ABIS" )
dev.off()

rownames(decon.df)
decon.df <- decon.df[c(1,3,4,5,7,8,10,11,13,15,17:23,29), ]
decon.df <- rbind(decon.df, colSums(decon.df[c("Bnsm", "Bsm"),]) )
decon.df <- rbind(decon.df, colSums(decon.df[c("CD8Tcm", "CD8Tem"),]) )
decon.df <- rbind(decon.df, colSums(decon.df[c("Mono_I", "Mono_NC"),]) )
rownames(decon.df)[19:21] <- c("Bm", "CD8Tm", "Mono_NC")
decon.df <- decon.df[c(1,2,5,8:11,14:21),]
rownames(decon.df)[12]  <- "Tgd"

pal.celltype <- brewer.pal(8, "Accent")
pal.celltype <- colorRampPalette(pal.celltype)(15)  
colnames(decon.df) <- pdata.benchmark$id
mat_col <- data.frame(group = pdata.benchmark$cell)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata.benchmark$cell)
pheatmap(decon.df[order(rownames(decon.df)),order(pdata.benchmark$cell) ],
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = NA, 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "ABIS" )
dev.off()







### SVM for topvargenes signatures####
for (i in c(100, 250, 500, 1000, 2000, 3000)) {
  keep <- intersect(rownames(sig.mat[topVarGenes[1:i],]), rownames(all.tpm))
  tpm.decon <- all.tpm[keep,]
  cibersort.df  <- cibersort(sig.mat[topVarGenes[1:i],], as.matrix(tpm.decon))
  cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))
  
  pal.celltype <- brewer.pal(8, "Set1")
  pal.celltype <- colorRampPalette(pal.celltype)(29)  
  colnames(cibersort.df) <- all.pdata$id
  mat_col <- data.frame(group = all.pdata$cell)
  rownames(mat_col) <- colnames(cibersort.df)
  mat_colors <- list(group = pal.celltype)
  names(mat_colors$group) <- unique(all.pdata$cell)
  #cibersort.df <- round(cibersort.df*100)
  pheatmap(cibersort.df[order(rownames(cibersort.df)),order(all.pdata$cell)], 
           cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
           show_colnames     = F,show_rownames     = T, border_color = "black", 
           annotation_col    = mat_col,annotation_colors = mat_colors,
           main = paste("RF meanGini", i , sep = " ") )
  print(paste(i, "Iteration Done", sep = " "))
}





#PCA
library(limma)
plotMDS(log2(sig.mat[rownames(meanGini)[1:500],]+1),  gene.selection="common", labels = colnames(sig.mat),
        col=pal.celltype[colnames(sig.mat)])
legend("top", legend=levels(factor(all.pdata$cell)), text.col=pal.celltype,
       bg="white", cex=0.7)



###in silico ####
library(reshape2)
library(ggplot2)
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
niteration <- 2280

all.pred.df <- data.frame(Pred = c(), Act = c(), CellType = c())
for (i in unique(all.pdata$cell)) {
  pred.act.df <- c()
  print(i)
  for (j in c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1)) {
    idx <- grep(i, all.pdata$cell)
    idx2 <- grep(i, all.pdata$cell, invert = T)
    simulation <- all.tpm[,idx2] * (1-j) + rowMeans( all.tpm[,idx]) * j #simulate the matrix
    simu_pred <- c()
    
    #SVM
    keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(simulation))
    simulation <- simulation[keep,]
    simu_pred  <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(simulation[,] , QN = F))
    #simu_pred  <- melt(simu_pred[,1:(ncol(simu_pred)-3)])
    simu_pred[,i]
    
    print(j); print("Done")
    tmp <- data.frame(Pred = simu_pred[,i], Act = rep(j, nrow(simu_pred)) , stringsAsFactors = F)
    pred.act.df <- rbind(pred.act.df, tmp)
  }
  pred.act.df$CellType <- rep(i, nrow(pred.act.df))
  all.pred.df <- rbind(all.pred.df, pred.act.df)
} 

#png("./plots/13Boxplot_PredAct_PMN.png", width = 4, height = 4, res = 300, units = "in")
#boxplot(pred.act.df$Pred ~ factor(pred.act.df$Act ,levels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1)),
#        main = "PMN", xlab = "Actual", ylab = "Prediction", col = "deepskyblue")
#dev.off()

load(file = "./dataset/COVID/insilico_Prediction_df.rds")
all.pred.df$Act <- factor(all.pred.df$Act ,  levels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1))
all.pred.df$CellType <- gsub("Mono_C", "cMono", all.pred.df$CellType)
all.pred.df$CellType <- gsub("Mono_I", "iMono", all.pred.df$CellType)
all.pred.df$CellType <- gsub("Mono_NC", "ncMono", all.pred.df$CellType)
all.pred.df$CellType <- gsub("VD2", "Tgd", all.pred.df$CellType)

all.pred.df$CellType <- factor(all.pred.df$CellType, levels = c("Erythroblast", "Megakaryocyte", "Basophils", 
                                                                "matureNeutrophil", "PMN", 
                                                                "cMono", "iMono", "ncMono", "mDC", "pDC",
                                                                "Bnaive", "Bsm", "Bnsm", "Bex", "Plasmablasts",
                                                                "CD4Tnaive", "CD4Teff", "Treg", "Th1", "Th2", "Th17", "Tfh", 
                                                                "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem","NK", "Tgd", "MAIT")
)



all.pred.df%>%
  ggplot( aes(x=Act, y=Pred)) +
  geom_boxplot(outlier.size=-1, aes( fill = "deepskyblue")) + ylim(0,1)+
  #geom_jitter(color="black", size=1, alpha=0.2) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("in silico Simulation ") + scale_fill_manual(values=c("deepskyblue"))+
  facet_wrap("CellType")+
  ylab("Percentage")+ 
  geom_abline(intercept = -0.1, slope = 0.1, linetype = 3,color = "darkblue")+theme_classic()

save(all.pred.df, file = "./dataset/COVID/insilico_Prediction_df.rds")





###20241113gene heatmap #####
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/gene_length.rds")
load(file = "./dataset/COVID/RF_meanGini.rds")
#load(file = "./dataset/COVID/RF_meanGini2.rds")
colnames(sig.mat)[c(17,18,19,29)] <- c("cMono", "iMono", "ncMono", "Tgd")
niteration <- 1220
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(29)  
names(pal.celltype) <-  c("Erythroblast", "Megakaryocyte", "Basophils", 
                          "matureNeutrophil", "PMN", 
                          "cMono", "iMono", "ncMono", "mDC", "pDC",
                          "Bnaive", "Bsm", "Bnsm", "Bex", "Plasmablasts",
                          "CD4Tnaive", "CD4Teff", "Treg", "Th1", "Th2", "Th17", "Tfh", 
                          "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem","NK", "Tgd", "MAIT")


sig.mat <- sig.mat[rownames(meanGini)[1:niteration ],]
mat_col <- data.frame(group = colnames(sig.mat))
rownames(mat_col) <- colnames(sig.mat)
mat_colors <- list(group = pal.celltype)
#names(mat_colors$group) <- colnames(sig.mat)
breaksList = seq(-1.5, 1.5, by = 0.1)
pheatmap(log10(sig.mat+1), 
         cluster_cols = F,cluster_rows = T, display_numbers = F, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList,
         show_colnames     = T,show_rownames     = F, border_color = NA, 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "1220 signature" )


factor(colnames(sig.mat), levels = c("Erythroblast", "Megakaryocyte", "Basophils", 
                                     "matureNeutrophil", "PMN", 
                                     "cMono", "iMono", "ncMono", "mDC", "pDC",
                                     "Bnaive", "Bsm", "Bnsm", "Bex", "Plasmablasts",
                                     "CD4Tnaive", "CD4Teff", "Treg", "Th1", "Th2", "Th17", "Tfh", 
                                     "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem","NK", "Tgd", "MAIT"))
cell.order <- c("Erythroblast", "Megakaryocyte", "Basophils", 
                "matureNeutrophil", "PMN", 
                "cMono", "iMono", "ncMono", "mDC", "pDC",
                "Bnaive", "Bsm", "Bnsm", "Bex", "Plasmablasts",
                "CD4Tnaive", "CD4Teff", "Treg", "Th1", "Th2", "Th17", "Tfh", 
                "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem","NK", "Tgd", "MAIT")
sig.mat <- sig.mat[,match(cell.order,  colnames(sig.mat))]

pal.bloodcell <- brewer.pal(3, "Accent")
mat_col <- data.frame(group = c(rep("Myeloid", 10), rep("B Lymphoid", 5), rep("T Lymphoid", 14)))
rownames(mat_col) <- colnames(sig.mat)
mat_colors <- list(group = c(rep("#7FC97F", 10), rep("#BEAED4", 5), rep("#FDC086", 14)))
names(mat_colors$group) <- mat_col$group
breaksList = seq(-1.5, 1.5, by = 0.1)
pheatmap(log10(sig.mat+1), 
         cluster_cols = F,cluster_rows = T, display_numbers = F, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList,
         show_colnames     = T,show_rownames     = F, border_color = NA, 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Blood Cell Type Signature" )
dev.off()



###Overlap of MELODY, ABIS, CIBERSORT####
sig.mat
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
sig_matrix <- read.delim(sig_matrix)
cibersort.gene <- sig_matrix$Gene.symbol

melody.gene <- rownames(meanGini)[1:niteration]
melody.gene <- gsub("_", "-", melody.gene)
intersect(melody.gene, cibersort.gene)

abis.gene <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt")
abis.gene <- rownames(abis.gene)
intersect(intersect(melody.gene, abis.gene), cibersort.gene)
intersect(melody.gene, abis.gene)
intersect(cibersort.gene, abis.gene)


#venn plot 

library(VennDiagram)
venn.diagram(
  x = list(cibersort.gene, abis.gene, melody.gene),
  category.names = c("CIBERSORT" , "ABIS" , "MELODY"),
  filename = "./dataset/COVID/13Venn_3algorithms.png",
  output=F
)




###SVM signatures####
keep <- intersect(rownames(sig.mat[sig.feature,]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[sig.feature,], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))

pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(29)  
colnames(cibersort.df) <- all.pdata$id
mat_col <- data.frame(group = all.pdata$cell)
rownames(mat_col) <- colnames(cibersort.df)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(all.pdata$cell)
#cibersort.df <- round(cibersort.df*100)
pheatmap(cibersort.df[order(rownames(cibersort.df)),order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = paste("RF meanGini2", "3000" , sep = " ") )





###1a. Overlap signatures####
ov <- intersect(genelist, rownames(abis_sig))
sig.mat <- sig.mat[ov,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[33] <- "VD2_925L2"
pdata$V5[89] <- "VD2_G4YW2"
pdata$V5[61] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Overlap signature" )



###1b. Merge signatures####
library(nnls)
library(pheatmap)
library(RColorBrewer)
abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./PBMC_GEO.RData")
load(file = "./dataset/COVID/PBMC_GEO_signature.RData")

ov <- unique(c(genelist, rownames(abis_sig)))
sig.mat <- sig.mat[ov,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[33] <- "VD2_925L2"
pdata$V5[89] <- "VD2_G4YW2"
pdata$V5[61] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Merge signature" )



###1d. RandomForest signatures####
library(nnls)
library(pheatmap)
library(RColorBrewer)
#abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/all.data.RData")

tpm <- all.tpm
#sig.mat <- sig.mat[RF.genelist,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[33] <- "VD2_925L2"
pdata$V5[89] <- "VD2_G4YW2"
pdata$V5[61] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "RandomForest signature" )



###1e. cutoff signatures####
library(nnls)
library(pheatmap)
library(RColorBrewer)
abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./PBMC_GEO.RData")
load(file = "./dataset/COVID/PBMC_GEO_signature.RData")

sig.mat <- sig.mat[cutoff.genelist,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[33] <- "VD2_925L2"
pdata$V5[89] <- "VD2_G4YW2"
pdata$V5[61] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Cutoff signature" )


library(ggplot2)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

auc <- data.frame(cell = c(), perc = c(), sd = c())
for (i in colnames(pdata)[5:ncol(pdata)]) {
  df2 <- data_summary(pdata, varname= i, 
                      groupnames=c("cell"))
  colnames(df2)[2] <- "perc"
  auc <- rbind(auc, df2[grep(i, df2$cell),])
}
ggplot(auc, aes(x=cell, y=perc, fill=cell)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=perc-sd, ymax=perc+sd), width=.2,
                position=position_dodge(.9)) +
  labs(title="Cutoff signature", x="Cells", y = "Perc.")+
  theme_classic() +
  scale_fill_manual(values=c(pal.celltype))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")



###1e. genelist signatures####
library(nnls)
library(pheatmap)
library(RColorBrewer)
abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./PBMC_GEO.RData")
load(file = "./dataset/COVID/PBMC_GEO_signature.RData")

sig.mat <- sig.mat[genelist,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[33] <- "VD2_925L2"
pdata$V5[89] <- "VD2_G4YW2"
pdata$V5[61] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "genelist signature" )




###1f. abis signatures####
library(nnls)
library(pheatmap)
library(RColorBrewer)
abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./PBMC_GEO.RData")
load(file = "./dataset/COVID/PBMC_GEO_signature.RData")

sig.mat <- sig.mat[rownames(abis_sig),]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[33] <- "VD2_925L2"
pdata$V5[89] <- "VD2_G4YW2"
pdata$V5[61] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "abis signature" )


library(ggplot2)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

auc <- data.frame(cell = c(), perc = c(), sd = c())
for (i in colnames(pdata)[5:ncol(pdata)]) {
  df2 <- data_summary(pdata, varname= i, 
                      groupnames=c("cell"))
  colnames(df2)[2] <- "perc"
  auc <- rbind(auc, df2[grep(i, df2$cell),])
}
ggplot(auc, aes(x=cell, y=perc, fill=cell)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=perc-sd, ymax=perc+sd), width=.2,
                position=position_dodge(.9)) +
  labs(title="ABIS signature", x="Cells", y = "Perc.")+
  theme_classic() +
  scale_fill_manual(values=c(pal.celltype))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")




###1g. RFovCutoff signatures####
library(nnls)
library(pheatmap)
library(RColorBrewer)
abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./PBMC_GEO.RData")
load(file = "./dataset/COVID/PBMC_GEO_signature.RData")

sig.mat <- sig.mat[intersect(cutoff.genelist, RF.genelist),]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[33] <- "VD2_925L2"
pdata$V5[89] <- "VD2_G4YW2"
pdata$V5[61] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "RFovCutoff signature" )




###2. RandomForest  signatures####
library(nnls)
library(pheatmap)
library(RColorBrewer)
#load("./dataset/COVID/PBMC_GEO_sigMat.rds")
#load(file = "./PBMC_GEO.RData")
#load(file = "./dataset/COVID/RF_sig.RData")
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/all.data.RData")

tpm <- all.tpm
#sig.mat <- sig.mat[RF.genelist157,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

pdata <- all.pdata
pdata <- cbind(pdata, as.data.frame(t(decon.all.df)))


library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
idx <-  grep("PBMC|Th1Th", pdata$cell, invert = T)
pdata <- pdata[idx,]
decon.all.df <- decon.all.df[,idx]

decon.heatmap <- round(decon.all.df,2)
pdata$V5 <- paste(pdata$cell, pdata$sample, sep = "_")
pdata$V5[7] <- "VD2_DZQV2"
pdata$V5[32] <- "VD2_925L2"
pdata$V5[86] <- "VD2_G4YW2"
pdata$V5[59] <- "VD2_9JD42"
colnames(decon.heatmap) <- pdata$V5
mat_col <- data.frame(group = pdata$cell)
rownames(mat_col) <- colnames(decon.heatmap)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdata$cell)
pheatmap(decon.heatmap[order(rownames(decon.heatmap)),order(pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "RandomForest 157 signature" )


library(ggplot2)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

auc <- data.frame(cell = c(), perc = c(), sd = c())
for (i in colnames(pdata)[5:ncol(pdata)]) {
  df2 <- data_summary(pdata, varname= i, 
                      groupnames=c("cell"))
  colnames(df2)[2] <- "perc"
  auc <- rbind(auc, df2[grep(i, df2$cell),])
}
ggplot(auc, aes(x=cell, y=perc, fill=cell)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=perc-sd, ymax=perc+sd), width=.2,
                position=position_dodge(.9)) +
  labs(title="RandomForest 500 signature", x="Cells", y = "Perc.")+
  theme_classic() +
  scale_fill_manual(values=c(pal.celltype))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")



###3. Markers of RF####
library(nnls)
library(pheatmap)
library(RColorBrewer)
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./PBMC_GEO.RData")
load(file = "./dataset/COVID/RF_sig.RData")

sig.mat <- sig.mat[RF.genelist157,]
sig.mat <- sig.mat[complete.cases(sig.mat),]
sig.mat[,grep("PBMC", colnames(sig.mat))]

mat_col <- data.frame(group = colnames(sig.mat))
rownames(mat_col) <- colnames(sig.mat)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- colnames(sig.mat)
breaksList = seq(-2, 2, by = 0.1)
pheatmap(log10(sig.mat+1), 
         cluster_cols = F,cluster_rows = T, display_numbers = F, scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList,
         show_colnames     = T,show_rownames     = T, border_color = NA, 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "RandomForest 157 signature" )

###3. in silico titration####
library(nnls)
library(pheatmap)
library(RColorBrewer)
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./PBMC_GEO.RData")
load(file = "./dataset/COVID/RF_sig.RData")

sig.mat <- sig.mat[RF.genelist157,]
sig.mat <- sig.mat[complete.cases(sig.mat),]

cell.name <- c()
for (i in 1:nrow(sig.mat)) {
  cell.name <- c(cell.name,names(which.max(sig.mat[i,])))
}
cell.manifest <- data.frame(gene = rownames(sig.mat), cellname = cell.name)
table(cell.manifest$cellname)






###5. validation cohort####
load( file = "./dataset/COVID/PBMC_GEO_validation.RData")
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/RF_meanGini.rds")
tpm
pdata
niteration <- 1220
rownames(tpm) <- gsub("-", "_", rownames(tpm))
keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(tpm))
tpm.decon <- tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))

pdata <- cbind(pdata, t(cibersort.df))
