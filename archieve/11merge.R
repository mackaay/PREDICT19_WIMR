library(edgeR)
library(matrixStats)
library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(33) 
load(file = "./dataset/COVID/gene_length.rds")


###1. merge datasets####
load(file = "./dataset/COVID/PBMC_GEO_ABIS.RData")
idx <- grep("PBMC|Th1Th17", pdata$cell, invert = T)
all.tpm <- tpm[,idx] 
all.pdata <- pdata[idx,]
all.count <- counts[,idx] 



load(file = "./dataset/COVID/PBMC_GEO_GMDSC_GSE178824.RData")
pdata
idx <-  grep("0", pdata$V4, invert = T)
all.tpm <- cbind(all.tpm, tpm[,idx])
all.count <- cbind(all.count, counts[,idx])
colnames(pdata) <- c("id", "cell", "sample", "V5")
all.pdata <- rbind(all.pdata, pdata[idx,])


load(file = "./dataset/COVID/PBMC_GEO_MMDSC_GSE199286.RData")
pdata
idx <-  grep("neg", pdata$V4)
all.tpm <- cbind(all.tpm, tpm[,idx])
all.count <- cbind(all.count, counts[,idx])
colnames(pdata) <- c("cell", "id", "sample", "V5")
all.pdata <- rbind(all.pdata, pdata[idx,])
all.pdata$cell[111:115] <- "M_MDSC"


load(file = "./dataset/COVID/PBMC_GEO_monoMMDSC_GSE183854.RData")
pdata
all.tpm <- cbind(all.tpm, tpm[,])
all.count <- cbind(all.count, counts[,])
colnames(pdata) <- c("V5", "id", "sample", "cell")
all.pdata <- rbind(all.pdata, pdata[,])


load(file = "./dataset/COVID/PBMC_GEO_neutrophil_GSE163533.RData")
pdata
all.tpm <- cbind(all.tpm, tpm[,])
all.count <- cbind(all.count, counts[,])
colnames(pdata) <- c("V5", "id", "sample", "cell")
all.pdata <- rbind(all.pdata, pdata[,])


load( file = "./dataset/COVID/PBMC_GEO_neutrophil_GSE66895.RData")
pdata
all.tpm <- cbind(all.tpm, tpm[,])
all.count <- cbind(all.count, counts[,])
colnames(pdata) <- c("V5", "id", "sample", "cell")
pdata$cell <- "neutrophil"
all.pdata <- rbind(all.pdata, pdata[,])


load(file = "./dataset/COVID/PBMC_GEO_PMN_GSE163834.RData")
pdata
idx <-  grep("healthy", pdata$V7)
all.tpm <- cbind(all.tpm, tpm[,idx])
all.count <- cbind(all.count, counts[,idx])
colnames(pdata) <- c("V5", "id", "cell", "sample")
all.pdata <- rbind(all.pdata, pdata[idx,])






####removing cells #####
#remove monocyte
idx <- grep("HD[1-3]", all.pdata$sample, invert = T)
all.tpm <- all.tpm[,idx]
all.count <- all.count[,idx]
all.pdata <- all.pdata[idx,]

# short reads 
idx <- grep("Peripheral", all.pdata$sample, invert = T)
all.tpm <- all.tpm[,idx]
all.count <- all.count[,idx]
all.pdata <- all.pdata[idx,]

#neutrophil
idx <- grep("JB[6-8]", all.pdata$sample, invert = T)
all.tpm <- all.tpm[,idx]
all.count <- all.count[,idx]
all.pdata <- all.pdata[idx,]


#
idx <- grep("healthy|KR|MS", all.pdata$sample, invert = T)
all.tpm <- all.tpm[,idx]
all.count <- all.count[,idx]
all.pdata <- all.pdata[idx,]



if (F){
  idx <- grep("Progenitor", all.pdata$cell, invert = T)
  all.tpm <- all.tpm[,idx]
  all.pdata <- all.pdata[idx,]
  all.count <- all.count[,idx]
  
  #M-MDSC
  idx <- grep("GSM5968715|GSM5968716", all.pdata$id, invert = T)
  all.tpm <- all.tpm[,idx]
  all.count <- all.count[,idx]
  all.pdata <- all.pdata[idx,]
  idx <- grep("M_MDSC", all.pdata$cell, invert = T)
  all.tpm <- all.tpm[,idx]
  all.count <- all.count[,idx]
  all.pdata <- all.pdata[idx,]
  
  #remove G_MDSC 
  idx <- grep("G_MDSC", all.pdata$cell, invert = T)
  all.tpm <- all.tpm[,idx]
  all.count <- all.count[,idx]
  all.pdata <- all.pdata[idx,]

  #LD Neutrophil
  idx <- grep("Neutrophils", all.pdata$cell, invert = T)
  all.tpm <- all.tpm[,idx]
  all.count <- all.count[,idx]
  all.pdata <- all.pdata[idx,]
}



table(all.pdata$cell)
all.pdata$cell <- gsub("Neutrophils|bandNeutrophil", "PMN", all.pdata$cell)
#all.pdata$cell <- gsub("neutrophils", "neutrophil", all.pdata$cell)



####TMM normalization####
library(edgeR)
y <- DGEList(counts = all.count[,], samples = all.pdata,
             group  = all.pdata$cell)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 4)

tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
all.tpm <- round(tpm, 4)


####Adding Erythro####
load(file = "./dataset/COVID/Erythro_Mega_Neutro.RData")
keep <- intersect(rownames(all.tpm), rownames(Erythro_Mega_Neutro))
all.tpm <- all.tpm[keep,]
Erythro_Mega_Neutro <- Erythro_Mega_Neutro[keep,]
all.tpm <- cbind(all.tpm, Erythro_Mega_Neutro)

all.pdata
pdata$sample <- pdata$id
pdata$V5 <- NA
all.pdata <- rbind(all.pdata, pdata)


#remove some erythro and band neutrophil , outlier
idx <- grep("band", all.pdata$cell, invert = T)
all.pdata <- all.pdata[idx, ]
all.tpm <- all.tpm[,idx]
idx <- grep("erythroblast.1|erythroblast.2|erythroblast.3|erythroblast.4", all.pdata$sample, invert = T)
all.pdata <- all.pdata[idx, ]
all.tpm <- all.tpm[,idx]


####changing name PMN####
#all.pdata$cell <- gsub("Neutrophils|bandNeutrophil", "PMN", all.pdata$cell)
#all.pdata$cell <- gsub("G.MDSCs", "PMN", all.pdata$cell)


####batch correction ######
#devtools::install_github("zhangyuqing/sva-devel")
library(sva)
all.pdata$batch <- 1
all.pdata$batch[grep("GSM", all.pdata$id, invert = T)] <- 2

adjusted.tpm <- ComBat_seq(all.tpm, batch = all.pdata$batch)




###Renaming the cell types #########
table(all.pdata$cell)
all.pdata$cell <- gsub("Mono_C", "cMono", all.pdata$cell)
all.pdata$cell <- gsub("Mono_I", "iMono", all.pdata$cell)
all.pdata$cell <- gsub("Mono_NC", "ncMono", all.pdata$cell)
all.pdata$cell <- gsub("matureNeutrophil", "Neutrophil", all.pdata$cell)
all.pdata$cell <- gsub("VD2", "Tgd", all.pdata$cell)
all.pdata$cell <- gsub("M_MDSC", "MMDSC", all.pdata$cell)
all.pdata$cell <- gsub("PMN", "PMNMDSC", all.pdata$cell)
all.pdata$cell <- factor(all.pdata$cell, levels = c("Progenitor", "Erythroblast", "Megakaryocyte", "Basophils", 
                                                    "Neutrophil", "PMNMDSC", 
                                                    "cMono", "iMono", "ncMono", "MMDSC","mDC", "pDC",
                                                    "Bnaive", "Bsm", "Bnsm", "Bex", "Plasmablasts",
                                                    "CD4Tnaive", "CD4Teff", "Treg", "Th1", "Th2", "Th17", "Tfh", 
                                                    "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem","NK", "Tgd", "MAIT"))

library(ggplot2)
pal.celltype <- colorRampPalette(brewer.pal(8, "Set1"))(31) 
names(pal.celltype) <- levels(all.pdata$cell)
save(pal.celltype, file = "./dataset/COVID/pal.celltype.rds")

###2. mito RPL genes ####
idx <- grep("^MT-", rownames(adjusted.tpm), value = T)
all.pdata$mito_perc <-colSums(adjusted.tpm[idx,])/colSums(adjusted.tpm)

idx <- grep("^RPL[0-9]", rownames(adjusted.tpm), value = T)
all.pdata$ribo_perc <-colSums(adjusted.tpm[idx,])/colSums(adjusted.tpm)

ggplot(all.pdata, aes(x=cell, y=mito_perc, color=cell)) +
  geom_boxplot() +scale_color_manual(values = pal.celltype) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  theme(legend.position = "none")
ggplot(all.pdata, aes(x=cell, y=ribo_perc, color=cell)) +
  geom_boxplot()+scale_color_manual(values = pal.celltype)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  theme(legend.position = "none")+
  labs(x = "Cell Type",                                              # labelling x axis
       y = "Ribosomal RNA Percentage",                                        # labeling y axis
       title = "",        # title
       fill = "") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme( #axis.text.x=element_blank(),
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )

pdf("./plots/11Boxplot_riboRNA.pdf", width = 6, height = 4)
ggplot(all.pdata, aes(x=cell, y=ribo_perc, color=cell)) +
  geom_boxplot()+scale_color_manual(values = pal.celltype)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  theme(legend.position = "left")+
  labs(x = "Cell Type",                                              # labelling x axis
       y = "Ribosomal RNA Percentage",                                        # labeling y axis
       title = "",        # title
       fill = "") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme( #axis.text.x=element_blank(),
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  )
dev.off()

#######20250116publication##########
#PCA
plotMDS(log2(all.tpm+1), top=1000, gene.selection="common", labels = all.pdata$cell,
                col=pal.celltype[factor(all.pdata$cell)], pch = 19)
legend("top", legend=levels(factor(all.pdata$cell)), text.col=pal.celltype,
       bg="white", cex=0.5)
plotMDS(log2(all.tpm+1), top=1000, gene.selection="common", 
        col=pal.celltype[factor(all.pdata$cell)], pch = 19)
plotMDS(log2(adjusted.tpm+1), top=2000, gene.selection="common", labels = all.pdata$cell,
        col=pal.celltype[factor(all.pdata$cell)])

plotMDS(log2(adjusted.tpm+1), top=2000, gene.selection="common", labels = all.pdata$sample,
        col=pal.celltype[factor(all.pdata$sample)])
legend("top", legend=levels(factor(all.pdata$cell)), text.col=pal.celltype,
       bg="white", cex=0.7)
dev.off()


pdf("./plots/11PCA_celltype.pdf", width = 4, height = 4)
plotMDS(log2(adjusted.tpm+1), top=2000, gene.selection="common", 
        col=pal.celltype[factor(all.pdata$cell)], pch = 19)
dev.off()


#sample distance cluster 
cell_dist <- dist(t(log2(all.tpm[grep("MIR[0-9]|^MT-|^RPL[0-9]", rownames(all.tpm), invert = T),]+1)))
cell_hclust <- hclust(cell_dist, method = "ward.D")
# The default `plot()` function can be used to produce a simple dendrogram
plot(cell_hclust, label = all.pdata$cell, type = "triangle", ylab = "Height")

cell_dist <- dist(t(log2(adjusted.tpm[grep("MIR[0-9]|^MT-|^RPL[0-9]", rownames(adjusted.tpm), invert = T),]+1)))
cell_hclust <- hclust(cell_dist, method = "ward.D")
# The default `plot()` function can be used to produce a simple dendrogram
plot(cell_hclust, label = all.pdata$cell, type = "triangle", ylab = "Height")

#nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),   cex = 0.7, col = "blue")
#plot(cell_hclust, label = all.pdata$cell, hang = -1, cex = 0.6,nodePar = nodePar, horiz = TRUE)
#plot(cell_hclust,type = "triangle", ylab = "Height")

library("ape")
cell_hclust$labels <- as.character( all.pdata$cell)
colors = brewer.pal(8, "Set2")[4:8]
clus4 = cutree(cell_hclust, 5)
plot(as.phylo(cell_hclust), type = "fan", cex = 1.2, tip.color = colors[clus4],
     no.margin = T)
dev.off()

pdf("./plots/11Phylo_celltype.pdf", width = 9, height = 9)
plot(as.phylo(cell_hclust), type = "fan", cex = 1, tip.color = colors[clus4],
     no.margin = T)
dev.off()



rownames(adjusted.tpm)
adjusted.tpm <- adjusted.tpm[grep("\\.", rownames(adjusted.tpm), invert = T),]
rownames(adjusted.tpm) <- gsub("-", "_", rownames(adjusted.tpm))
#all.tpm <- all.tpm[grep("RPL[0-9]", rownames(all.tpm), invert = T),]


save(adjusted.tpm, all.pdata, file ="./dataset/COVID/all_adjusted.RData")





###3. PMN vs neutrophil, M-MDSC vs mono####
cellType <- factor(all.pdata$cell)
design <- model.matrix(~0+cellType, data=all.pdata)
colnames(design) <- c(levels(cellType))
contMatrix <- makeContrasts(PMNMDSC-Neutrophil,
                            levels=design)
contMatrix

v <- voom(adjusted.tpm, design, plot=TRUE)
v
vfit <- lmFit(log2(all.tpm[grep("MIR[0-9]|^MT_|^RPL[0-9]", rownames(all.tpm), invert = T),]+1), design)
vfit <- contrasts.fit(vfit, contrasts=contMatrix)
efit <- eBayes(vfit)
summary(decideTests(efit))
DEG <- topTable(efit, num=Inf, coef=1)
head(DEG)

rownames(DEG)[which(abs(DEG$logFC) > 4 & DEG$adj.P.Val < 0.0001)]
which(abs(DEG$logFC) > 4 & DEG$adj.P.Val < 0.0001)
DEG[c("HDC", "CD33"),]
all.tpm[c("HDC", "CD33"),grep("PMN|G_MDSCs|neutro", all.pdata$cell)]



#### top variable genes####
topVarGenes <- head(order(rowVars(as.matrix(adjusted.tpm)), decreasing = TRUE), 5000)
rownames(all.tpm)[topVarGenes]





matrix(1:20, nrow = 4)*c(2,2,3,3,0)
sweep(matrix(1:20, nrow = 4), MARGIN = 2, c(2,2,3,3,0), '*')
