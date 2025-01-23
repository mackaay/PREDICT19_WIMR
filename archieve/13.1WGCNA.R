library(dynamicTreeCut)
library(ape)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(RColorBrewer)

load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/gene_length.rds")
load(file = "./dataset/COVID/RF_meanGini.rds")
load(file = "./dataset/COVID/RF_meanGini2.rds")
load(file = "./dataset/COVID/all.data.RData")
niteration <- 1220
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(29)  



#####phylo and bootstrap#######
lcpm.cell <- c()
for (i in levels(as.factor(all.pdata$cell))  ) {
  print(i)
  sra.id <- all.pdata$id[grep(i, all.pdata$cell)]
  lcpm.cell <- cbind(lcpm.cell, rowMeans(all.tpm[, sra.id]) )
}
colnames(lcpm.cell) <- levels(as.factor(all.pdata$cell))
lcpm.cell[grep("CD74",rownames(lcpm.cell)) ,]
#ape distance and neighbour-joining 
dist.lcpm.cell <- dist.gene( t(lcpm.cell) , method = "pairwise") #Transpose
length(dist.lcpm.cell)
tre <- nj(dist.lcpm.cell)
class(tre) 
tre <- ladderize(tre)
tre # tells us what the tree will look like but doesn't show the actual construction
plot(tre, cex = 1.5)
title("A Simple NJ Tree")

#checking 
x <- as.vector(dist.lcpm.cell)
tre2 <- root(tre, out = 1)
tre2 <- ladderize(tre2)
y <- as.vector(as.dist(cophenetic(tre2)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree", 
     main="Is NJ appropriate?", pch=20, col="black", cex=3)
abline(lm(y~x), col="red")


myBoots <- boot.phylo(tre2, t(lcpm.cell), function(e) root(nj(dist.gene(e, method = "pairwise")),1))
myBoots
png("./plots/13.1Phylo_Cell.png", width = 100, height = 100, units = "mm", res = 300)
plot(tre2, show.tip=T, edge.width=2, cex = 1.2)
title("NJ tree + bootstrap values")
#nodelabels(myBoots, cex=.1)
dev.off()




#####1. modules of DEG#####
library(WGCNA)
enableWGCNAThreads(nThreads = 4)
#RNAseq_tumor = RNAseq_tumor[apply(RNAseq_tumor,1,function(x) sum(x==0))<ncol(RNAseq_tumor)*0.8,]
DEG <- rownames(meanGini)[1:1220]
keep <- intersect(DEG, rownames(all.tpm))
lcpm.DEG <- all.tpm[keep, ]
DEG.df <-  t(lcpm.DEG)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(DEG.df, powerVector = powers, verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='black')
abline(h=0.85,col='red')
sft$powerEstimate


#calculation of adjacency matrix
#beta = sft$powerEstimate
#a = s^beta
#dissimilarity measure
#w = 1-a
net = blockwiseModules(DEG.df, power = 9, maxBlockSize = 6000, 
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,saveTOMs = F,  verbose = 3)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
png("./plots/01.1Dengrogram_DEG.png", width = 10 ,height = 10,units = "cm",  res = 300)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
table(net$colors)
modulelabels <- net$colors
names(which(modulelabels == 0))
module.df <- data.frame(cluster = modulelabels, gene = names(modulelabels), color = mergedColors, 
                        stringsAsFactors = F)


##Complex Heatmap of GO from DEG Module 
library(clusterProfiler)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

GO.DEG <- c()
for (i in 0:11) {
  print(i)
  my.symbols <- module.df$gene[module.df$cluster == i]
  entrezid <- select(hs, 
                     keys = my.symbols,
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "SYMBOL")
  egoCC <- enrichGO(gene          = entrezid$ENTREZID,
                    OrgDb         = 'org.Hs.eg.db',
                    keyType = "ENTREZID",
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 5, maxGSSize = 1000,
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE, pool = T)
  egoCC <- head(egoCC)
  egoCC$clusterid <- rep(i, nrow(egoCC)  )
  print("GO done, merging...")
  GO.DEG <- rbind(GO.DEG, egoCC)
}
write.csv(GO.DEG, file = "./dataset/COVID/GO.DEG.csv" , col.names = T, row.names = F)
write.csv(module.df, file = "./dataset/COVID/DEG_module.csv", col.names = T, row.names = F)

GO.DEG[,c("Description", "clusterid")]
table(module.df$color, module.df$cluster)


#heatmap
module.df <- read.csv( file = "./dataset/COVID/DEG_module.csv")
library(ComplexHeatmap)
table(module.df$cluster)
module.df <- module.df[order(module.df$cluster, decreasing = F),] #Reordering the module cluster
lcpm.DEG <- lcpm.DEG[match(module.df$gene, rownames(lcpm.DEG)),]
all.pdata <- all.pdata[order(all.pdata$cell, decreasing = F),] # Reordering the cell type sample
lcpm.DEG <- lcpm.DEG[,match(all.pdata$id, colnames(lcpm.DEG))]
module.df$module <- paste0("module", module.df$cluster)
module.df$color <- as.character(module.df$color)
module.anno.col <- as.character(module.df$color)
names(module.anno.col) <- module.df$module

rowMax(lcpm.DEG)

DEG.label <- c("CCL5","KLRK1","IFNG","CD8B", "CD8A" ,"CXCR6", "CD4",
               "SELL", "CCR7", "LEF1", "CTLA4",
               "GZMH","FOXP3","GZMK","NCR3","KLRB1","TCL1A","CD22",
               "PDCD1","LAG3","TIGIT", "HAVCR2",  "EPHA1","IL1RL1" , 
               "CD79A", "CD79B", "CD19","MS4A1", #B cell
               "JCHAIN", "XBP1",  #plasma
               "CD1C", "CD1C", "CD1D", "CD1E" ,"CLEC9A", #DC
               "KNRD1", "XCL2", "XCL1", "NKG7", "NCAM1"  #NK
               )
eset <- log2(lcpm.DEG+1)
data.frame(all.pdata$id, colnames(eset))
all.pdata$cell <- factor(all.pdata$cell, levels = c("Erythroblast", "Megakaryocyte", "Basophils", 
                                     "matureNeutrophil", "PMN", 
                                     "Mono_C", "Mono_I", "Mono_NC", "mDC", "pDC",
                                     "Bnaive", "Bsm", "Bnsm", "Bex", "Plasmablasts",
                                     "CD4Tnaive", "CD4Teff", "Treg", "Th1", "Th2", "Th17", "Tfh", 
                                     "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem","NK", "VD2", "MAIT"))
eset <- eset[,order(all.pdata$cell)]
all.pdata <- all.pdata[order(all.pdata$cell),]
names(pal.celltype) <- levels(factor(all.pdata$cell))

circlize::colorRamp2(c(-2, 0, 12), c("blue", "white", "red"))
colorRampPalette(c("blue", "white", "red"))(50)
png("./plots/13.1Heatmap_ModuleDEG.png", width = 17, height = 27, res= 300, units = "cm")
Heatmap(eset, col = circlize::colorRamp2(c(0, 5, 15), c("blue", "white", "red")), 
        name = "Exp", row_split = module.df$module , column_km = 0, border = T,
        #row_split = rep(c("A", "B"), 9), column_split = rep(c("C", "D"), 12),
        cluster_rows = F, cluster_columns = F,
        row_title = "Module of Cell Type Signatures", row_title_side = "left",
        column_title = "Cell Types",
        show_column_names = FALSE, width = unit(8, "cm"),
        show_row_names = F,
        #top_annotation = HeatmapAnnotation(CellType = all.pdata$cell, col = list(CellType = pal.celltype)), 
        #left_annotation = rowAnnotation(Module = module.df$module, col = list(Module = module.anno.col)),
        heatmap_legend_param = list(title = "Log2(TPM)")) +
  rowAnnotation(link = anno_mark(at = which(rownames(lcpm.DEG) %in% DEG.label) , 
                                 labels = rownames(lcpm.DEG)[which(rownames(lcpm.DEG) %in% DEG.label) ], 
                                 labels_gp = gpar(fontsize = 7), 
                                 #padding = unit(1, "mm")
  ))
dev.off()
png("./plots/13.1Heatmap_ModuleDEG2.png", width = 17, height = 27, res= 300, units = "cm")
Heatmap(eset, col = circlize::colorRamp2(c(0, 5, 15), c("blue", "white", "red")), 
        name = "Exp", row_split = module.df$module , column_km = 0, border = T,
        #row_split = rep(c("A", "B"), 9), column_split = rep(c("C", "D"), 12),
        cluster_rows = F, cluster_columns = F,
        row_title = "Module of DEG", row_title_side = "left",
        column_title = "Cell Types",
        show_column_names = FALSE, width = unit(8, "cm"),
        show_row_names = F,
        top_annotation = HeatmapAnnotation(CellType = all.pdata$cell, col = list(CellType = pal.celltype )), 
        left_annotation = rowAnnotation(Module = module.df$module, col = list(Module = module.anno.col)),
        heatmap_legend_param = list(title = "Scaled expr")) 
dev.off()
table(sampleinfo$cell)



barplot(egoCC_sel, showCategory=10, title = 'GO CC of selected genes')
write.csv(egoCC_sel, file ='GOCC_sel.csv')



sampleTree = hclust(dist(lcpm.DEG), method = 'average')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
module.sample = cutreeDynamic(dendro = sampleTree,  deepSplit = 4, pamRespectsDendro = FALSE,
                              minClusterSize = 6)
#assign module colours
module.colours = labels2colors(modules)
#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='WGCNA of Prostate Cancer (TCGA)')
table(module.colours)




#####2. co-expressed module WGCNA ######


##DEG of PMN vs PMN-MDSC######
library(limma)
ltpm <- log2(all.tpm+1)

cellType <- factor(all.pdata$cell)
design <- model.matrix(~0+cellType, data=all.pdata)
colnames(design) <- c(levels(cellType))
fit <- lmFit(ltpm, design)
contMatrix <- makeContrasts(PMN-matureNeutrophil,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

DEG <- topTable(fit2, num=Inf, coef=1)
head(DEG)




grep("ARTC2|P2RX7|TBX21", rownames(sig.mat), value = T)
tmp <- cbind(all.pdata, t(all.tpm[c("P2RX7","TBX21"),]))
tmp %>%
  ggplot( aes(x=cell, y=P2RX7, fill=cell)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.4)

