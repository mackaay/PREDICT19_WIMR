library(edgeR)
library(matrixStats)
library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(32) 
load(file = "./dataset/COVID/gene_length.rds")
library(celldex)
ref <- BlueprintEncodeData()
ref

###Blueprint datasets ####
unique(colData(ref)$label.fine)
idx.tissue <- grep("CLP|CMP|GMP|HSC|MEP|MPP|mv |Chondrocytes|Mesangial|Preadipo|Erythr|Megaka|Eosino|Astro|Neuron|CD|Mono",colData(ref)$label.fine, invert = T)
idx.tissue <- grep("Neutrophil|Macrophage|Fibroblast|Melanocyte|Epithelial|^Endothelial",colData(ref)$label.fine)

target_tissue <- 2^assays(ref[,idx.tissue])[["logcounts"]] -1
target_tissue <- target_tissue[,order(colnames(target_tissue))]
colnames(target_tissue)
pdata <- colData(ref)[idx.tissue,]
keep <- intersect(rownames(pdata), colnames(target_tissue))
pdata <- pdata[keep, ]
target_tissue <- target_tissue[,keep]
pdataENCODE <- pdata
colnames(pdataENCODE) <- c("sample", "cell", "V5")
pdataENCODE$id <- rownames(pdataENCODE)
pdataENCODE <- as.data.frame(pdataENCODE)
#pdataENCODE$cell[grep("band", pdataENCODE$id)] <- "PMN"
pdataENCODE <- pdataENCODE[grep("band|bone.marrow", pdataENCODE$id, invert = T),]
target_tissue <- target_tissue[,grep("band|bone.marrow", pdataENCODE$id, invert = T)]

unique(pdataENCODE$cell)

###1. merge datasets####
load(file = "./dataset/COVID/PBMC_GEO_ABIS.RData")
idx <- grep("PBMC|Th1Th17", pdata$cell, invert = T)
all.tpm <- tpm[,idx] 
all.pdata <- pdata[idx,]
all.count <- counts[,idx] 
####removing cells #####
idx <- grep("Progenitor|Basophil|Mono", all.pdata$cell, invert = T)
all.tpm <- all.tpm[,idx]
all.pdata <- all.pdata[idx,]
all.count <- all.count[,idx]

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

all.pdata$cell <- gsub("Neutrophils", "PMN", all.pdata$cell)

####merging two datasets####
keep <- intersect(rownames(all.tpm), rownames(target_tissue))
all.tpm <- all.tpm[keep,]
target_tissue <- target_tissue[keep,]
all.tpm <- cbind(all.tpm, target_tissue)

all.pdata <- rbind(all.pdata, pdataENCODE)
all.pdata



#remove some erythro and band neutrophil
idx <- grep("mv ", all.pdata$cell, invert = T)
all.pdata <- all.pdata[idx, ]
all.tpm <- all.tpm[,idx]
idx <- grep("erythroblast.1|erythroblast.2|erythroblast.3|erythroblast.4", all.pdata$id, invert = T)
all.pdata <- all.pdata[idx, ]
all.tpm <- all.tpm[,idx]
idx <- grep("macrophage.1|macrophage.2|macrophage.3|macrophage.4", all.pdata$id, invert = T)
all.pdata <- all.pdata[idx, ]
all.tpm <- all.tpm[,idx]
idx <- grep("^dermis|glomerular|pericard", all.pdata$id, invert = T)
all.pdata <- all.pdata[idx, ]
all.tpm <- all.tpm[,idx]


####Remove batch effect ####
all.pdata$batch <- "ENCODE"
all.pdata$batch[grep("GSM", all.pdata$id)] <- "ABIS"
library(sva)
all.tpm <- ComBat_seq(all.tpm, batch = all.pdata$batch)


###2. PCA and Clusteringg#####
library(ggplot2)
unique(all.pdata$cell)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(34) 
pal.batch <- brewer.pal(3, "Dark2")


#tmp <- removeBatchEffect(all.tpm, all.pdata$batch)
#all.tpm <- tmp
#boxplot(log2(tmp+5))

#PCA
plotMDS(log2(all.tpm+1), top=1000, gene.selection="common", labels = all.pdata$cell,
        col=pal.celltype[factor(all.pdata$cell)])
plotMDS(tmp, top=1000, gene.selection="common", labels = all.pdata$cell,
        col=pal.batch[factor(all.pdata$batch)])
plotMDS(log2(target_tissue+1), top=1000, gene.selection="common", labels = pdataENCODE$cell,
        col=pal.celltype[factor(pdataENCODE$cell)])
#Cluster
cell_dist <- dist(t(log2(all.tpm[grep("MIR[0-9]|^MT-|^RPL[0-9]", rownames(all.tpm), invert = T),]+1)))
cell_hclust <- hclust(cell_dist, method = "ward.D")
# The default `plot()` function can be used to produce a simple dendrogram
plot(cell_hclust, label = all.pdata$cell, type = "triangle", ylab = "Height")

cell_dist <- dist(t(log2(target_tissue[grep("MIR[0-9]|^MT-|^RPL[0-9]", rownames(target_tissue), invert = T),]+1)))
cell_hclust <- hclust(cell_dist, method = "ward.D")
# The default `plot()` function can be used to produce a simple dendrogram
plot(cell_hclust, label = pdataENCODE$cell, type = "triangle", ylab = "Height")


library("ape")
cell_hclust$labels <- pdataENCODE$cell
colors = brewer.pal(5, "Set1")
clus4 = cutree(cell_hclust, 5)
plot(as.phylo(cell_hclust), type = "fan", cex = 0.7, tip.color = colors[clus4],
     no.margin = F)



###3. top variable genes####
library(randomForest)
library(caret)
library(edgeR)
library(matrixStats)
topVarGenes <- head(order(rowVars(as.matrix(all.tpm)), decreasing = TRUE), 5000)
rownames(all.tpm)[topVarGenes]


####First iteration ####
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
all.tpm <- all.tpm[grep("^MT_", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^RPL[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^MIR[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^LINC[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^RNU[0-9]*", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("\\.", rownames(all.tpm), invert = T),]
topVarGenes <- head(order(rowVars(as.matrix(all.tpm)), decreasing = TRUE), 5000)
rownames(all.tpm)[topVarGenes]
data <- as.data.frame(t(all.tpm[topVarGenes,]))
data$CellType <- all.pdata$cell
data$CellType <- as.factor(data$CellType)
table(data$CellType)

set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T) 
print(rf)
rf
varImpPlot(rf,
           sort = T,
           n.var = 30,
           main = "Top 30 - Variable Importance")
meanGini <- importance(rf)
meanGini <- meanGini[order(meanGini[,"MeanDecreaseAccuracy"], decreasing = T),]
meanGini2 <- importance(rf)[order(importance(rf)[,"MeanDecreaseGini"], decreasing = T),]

####sig.mat ######
sig.mat <- all.tpm[, ]
all.pdata <- cbind(all.pdata, t(sig.mat))
sig.mat <- aggregate(x = all.pdata[,c(6:ncol(all.pdata))],by = list(all.pdata$cell), FUN = median)
rownames(sig.mat) <- sig.mat$Group.1
sig.mat <- sig.mat[,-1]
sig.mat <- t(sig.mat)


####in silico #####
library(CIBERSORT)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(31) 

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:2000],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:2000],], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))

library(pheatmap)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(34)  
colnames(cibersort.df) <- all.pdata$id
mat_col <- data.frame(group = all.pdata$cell)
rownames(mat_col) <- colnames(cibersort.df)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(all.pdata$cell)
#cibersort.df <- round(cibersort.df*100)
pheatmap(cibersort.df[order(rownames(cibersort.df)),order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "RF 2000 CIBERSORT ENCODE" )

####Remove outliers samples####
tmp <- cibersort.df[,grep("Fibro|Endo|Epithe", all.pdata$cell)]
tmp.pdata <- all.pdata[grep("Fibro|Endo|Epithe", all.pdata$cell),]
pheatmap(tmp[order(rownames(tmp)),order(tmp.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "RF500 CIBERSORT ENCODE" )

tmp <- cibersort.df[,grep("Fibro", all.pdata$cell)]
rm.samples <- names(which(tmp["Fibroblasts",]<0.5))
tmp <- cibersort.df[,grep("Endothelial cells", all.pdata$cell)]
rm.samples <-  c(rm.samples, names(which(tmp["Endothelial cells",]<0.5)))
tmp <- cibersort.df[,grep("Epithelial cells", all.pdata$cell)]
rm.samples <-  c(rm.samples, names(which(tmp["Epithelial cells",]<0.5)))
tmp <- cibersort.df[,grep("Melanocytes", all.pdata$cell)]
rm.samples <-  c(rm.samples, names(which(tmp["Melanocytes",]<0.5)))
#tmp <- cibersort.df[,grep("Macrophage", all.pdata$cell)]
#rm.samples <-  c(rm.samples, names(which(tmp["Macrophages",]<0.5)))

####cleaning datasets####
idx <- which(! all.pdata$id %in% rm.samples)
all.pdata <- all.pdata[idx,]
all.tpm <- all.tpm[,idx]

topVarGenes <- head(order(rowVars(as.matrix(all.tpm)), decreasing = TRUE), 5000)
rownames(all.tpm)[topVarGenes]
data <- as.data.frame(t(all.tpm[topVarGenes,]))
data$CellType <- all.pdata$cell
data$CellType <- as.factor(data$CellType)
table(data$CellType)

set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T) 
print(rf)
rf
varImpPlot(rf,
           sort = T,
           n.var = 30,
           main = "Top 30 - Variable Importance")
meanGini <- importance(rf)
meanGini <- meanGini[order(meanGini[,"MeanDecreaseAccuracy"], decreasing = T),]
meanGini2 <- importance(rf)[order(importance(rf)[,"MeanDecreaseGini"], decreasing = T),]


####sig.mat ######
sig.mat <- all.tpm[, ]
all.pdata <- all.pdata[,1:5]
all.pdata <- cbind(all.pdata, t(sig.mat))
sig.mat <- aggregate(x = all.pdata[,c(6:ncol(all.pdata))],by = list(all.pdata$cell), FUN = median)
rownames(sig.mat) <- sig.mat$Group.1
sig.mat <- sig.mat[,-1]
sig.mat <- t(sig.mat)

table(colnames(sig.mat)[apply(sig.mat,1 ,which.max)])

####in silico #####
library(CIBERSORT)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(31) 

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:2000],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:2000],], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))

library(pheatmap)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(30)  
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
         main = "meanGini 2000 CIBERSORT ENCODE" )

keep <- intersect(rownames(sig.mat[rownames(meanGini2)[1:2000],]), rownames(all.tpm))
tpm.decon <- all.tpm[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini2)[1:2000],], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))

pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(30)  
colnames(cibersort.df) <- all.pdata$id
mat_col <- data.frame(group = all.pdata$cell)
rownames(mat_col) <- colnames(cibersort.df)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(all.pdata$cell)
#cibersort.df <- round(cibersort.df*100)
pheatmap(cibersort.df[order(rownames(cibersort.df)),order(all.pdata$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "meanGini2 2000 CIBERSORT ENCODE" )

save(meanGini, meanGini2, sig.mat, file = "./dataset/COVID/Tissue_sig.mat.RData")


####Iteration####
load(file = "./dataset/COVID/Tissue_sig.mat.RData")
accuracy <- data.frame(ngene =c(), Cor = c(), RMSE = c())
for (i in 2:200) {
  keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:(i*10)],]), rownames(all.tpm))
  tpm.decon <- all.tpm[keep,]
  cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:(i*10)],], as.matrix(tpm.decon))
  print((i*10))
  print(median(cibersort.df[,"Correlation"]))
  print(median(cibersort.df[,"RMSE"]) )
  tmp.acc <- data.frame(ngene = rep((i*10), nrow(cibersort.df)), 
                        Cor = cibersort.df[,"Correlation"], 
                        RMSE = cibersort.df[,"RMSE"])
  accuracy <- rbind(accuracy, tmp.acc)
  save(accuracy, file = "./dataset/COVID/SVM_accuracy_tissue.rds")
}

library(ggplot2)
library(dplyr)
accuracy$CellType <- rep(all.pdata$cell, nrow(accuracy)/nrow(all.pdata))
accuracy[nrow(accuracy),]
accuracy.df <- accuracy %>% 
  group_by(ngene) %>% 
  summarise(Correlation = median(Cor, na.rm = TRUE), RMSE= median(RMSE, na.rm = TRUE))
ggplot(data = accuracy.df, aes(x = ngene, y = Correlation))+
  geom_line()+ ylim(0.2,1)+ 
  geom_vline(xintercept=c(1700,2000), linetype="dotted", color = "darkred")
#geom_vline(xintercept=c(325),  color = "yellow", size = 25, alpha = 0.2)
ggplot(data = accuracy.df, aes(x = ngene, y = RMSE))+
  geom_line()+ ylim(0.25,1)+ 
  geom_vline(xintercept=c(1700,2000), linetype="dotted", color = "darkred")


###4. Melanoma PD1 ####
load(file = "./dataset/COVID/Tissue_sig.mat.RData")
counts <- read.delim("/datasets/work/hb-diab-cfdna/work/scratch/chenkai/melanoma_pd1/data/cancercell_normalized_counts_genenames.txt")
#counts <- read.delim("/datasets/work/hb-diab-cfdna/work/scratch/chenkai/melanoma_pd1/data/PD1-IPIPD1_counts.txt")
counts <- aggregate(x = counts[,4:ncol(counts)],by = list(counts$Gene), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
keep <- intersect(rownames(counts), names(gene.length))
counts <- counts[keep,]
gene.length <- gene.length[keep]
data.frame(names(gene.length), rownames(counts))

tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

pData.cancer <- read.csv("/datasets/work/hb-diab-cfdna/work/scratch/chenkai/melanoma_pd1/data/pData.csv", row.names = 1)
pData.cancer$Best.RECIST.response <- gsub(" ", "", pData.cancer$Best.RECIST.response )
pData.cancer$RNA.Sequencing <- gsub(" ", "", pData.cancer$RNA.Sequencing )
pData.cancer$Best.RECIST.response <- factor(pData.cancer$Best.RECIST.response, levels =  c("PD", "SD", "PR", "CR"))


decon.cibersort <- cibersort(sig.mat[rownames(meanGini2)[1:2000],grep("Plasma",colnames(sig.mat), invert = T)], as.matrix(tpm))
decon.cibersort <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])
keep <- intersect(colnames(decon.cibersort), rownames(pData.cancer))
decon.cibersort <- decon.cibersort[,keep]
pData.cancer <- pData.cancer[keep,]

mat_col <- data.frame(Response = pData.cancer$Best.RECIST.response, 
                      Time = pData.cancer$RNA.Sequencing)
rownames(mat_col) <- rownames(pData.cancer)
# List with colors for each annotation.
mat_colors <- list(Best.RECIST.response = brewer.pal(4, "Set2"),
                   Time = brewer.pal(3, "Dark2")[1:2])
names(mat_colors$Best.RECIST.response) <- unique(pData.cancer$Best.RECIST.response)
names(mat_colors$Time) <- unique(pData.cancer$RNA.Sequencing)
pheatmap(decon.cibersort[,order(pData.cancer$Best.RECIST.response)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Melanoma IPIPD1 cohort" )
decon.cibersort["Neutrophils",]


decon.cibersort_PRE <- decon.cibersort[,grep("PRE", pData.cancer$RNA.Sequencing)]
pData.cancer_PRE <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),]
mat_col <- data.frame(Response = pData.cancer_PRE$Best.RECIST.response)
rownames(mat_col) <- rownames(pData.cancer_PRE)
mat_colors <- list(Best.RECIST.response = brewer.pal(4, "Set2")
                   )
names(mat_colors$Best.RECIST.response) <- unique(pData.cancer_PRE$Best.RECIST.response)
pheatmap(decon.cibersort_PRE[,order(pData.cancer_PRE$Best.RECIST.response)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Melanoma IPIPD1 cohort" )





library(ggplot2)
library(ggpubr)
library(ggsignif)
data.frame(colnames(decon.cibersort), rownames(pData.cancer))
pData.cancer <- cbind(pData.cancer, t(decon.cibersort))
pData.cancer$Response <- "Response"
pData.cancer$Response[grep("PD|SD",pData.cancer$Best.RECIST.response)] <- "noResponse"

colnames(pData.cancer) <- gsub(" " , "_", colnames(pData.cancer))
p0 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Melanocytes, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Melanocytes PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.2)
p1 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=CD8Tem, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("CD8Tem PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.2)
p2 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Neutrophils, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Neutrophils PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.2)
p3 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Fibroblasts, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Fibroblasts PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.4)
p4 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Macrophages_M2, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Macrophage M2 PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.1)
p5 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Treg, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Treg PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.1)
ggarrange( p0, p1, p2, p3, p4,p5,
          labels = c("A", "B", "C", "D", "E","F"),
          ncol = 3, nrow = 2)


p0 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Melanocytes, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Melanocytes EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.1)
p1 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=CD8Tem, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("CD8Tem EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.3)
p2 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Neutrophils, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Neutrophils EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.1)
p3 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Fibroblasts, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Fibroblasts EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.4)
p4 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Macrophages_M2, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Macrophage M2 EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.1)
p5 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=Treg, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("Treg EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.1)
ggarrange( p0, p1, p2, p3, p4,p5,
           labels = c("A", "B", "C", "D", "E","F"),
           ncol = 3, nrow = 2)



###5. MIA scRNAseq  FFPE tissue####
#https://www.cell.com/cell-reports/fulltext/S2211-1247(24)00720-4
#remotes::install_github("omnideconv/immunedeconv")
library(ggplot2)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(readr)
library(tibble)
library(infercnv)
scale_to_million <- function(sample) {
  (sample / sum(sample)) * 1e6
}

ref.bulk <- log2(sig.mat[rownames(meanGini)[1:2000],]+1)

pdata.tmp <- data.frame(id = all.pdata$id, sample = all.pdata$sample, cell = all.pdata$cell)
counts.tmp <- as.matrix(log2(all.tpm[,]+1))
ref.bulk.all <- SummarizedExperiment(assays=list(logcounts=counts.tmp),
                                 colData=pdata.tmp)
ref.bulk.all
ref.bulk <- SummarizedExperiment(assays=list(logcounts=counts.tmp[rownames(meanGini)[1:2000],]),
                                 colData=pdata.tmp)
ref.bulk

pbmc.data <- read.csv("./dataset/COVID/scRNAseq_MIA/Melanoma_RNA_counts.csv", row.names = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "Patient6", min.cells = 3, min.features = 200)
pbmc
metadata <- read.csv("./dataset/COVID/scRNAseq_MIA/Melanoma_RNA_metadata.csv", row.names = 1)
keep <- intersect(rownames(metadata), colnames(pbmc))
metadata <- metadata[keep,]
pbmc <- pbmc[,keep]
data.frame(rownames(metadata), colnames(pbmc))
pbmc$MIA_ID <- metadata$mia_ID
pbmc$RECIST <- as.factor(metadata$RECIST)
pbmc$timepoint <- as.factor(metadata$timepoint)

pbmc.sce <- as.SingleCellExperiment(pbmc)
pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred.singler$labels)
pbmc$singler.fine <- pred.singler$labels
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=ref.bulk$cell )
table(pred.bulk$labels)
pbmc$bulk <- pred.bulk$labels
pred.bulk.all <- SingleR(test=pbmc.sce, ref=ref.bulk.all, labels=ref.bulk.all$cell )
table(pred.bulk.all$labels)
pbmc$bulk.all <- pred.bulk.all$labels
#pred.df <- data.frame(SingleR = pred.singler$labels, bulk = pred.bulk$labels,  freq = c(1))

pbmc <- pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() #features = rownames(pbmc.combined)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims = 30)
pbmc <- FindNeighbors(object = pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = c( 0.1, 0.3))
pbmc <- RunUMAP(pbmc, reduction = "pca",  dims = 1:20)
Idents(pbmc) <- pbmc[["RNA_snn_res.0.3"]] 
DimPlot(pbmc, label = TRUE, reduction = "umap") + NoLegend() +ggtitle("All")


pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(15) 
Idents(pbmc) <- 'RNA_snn_res.0.3'
p2 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(38) 
Idents(pbmc) <- 'singler.fine'
p3 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(29) 
Idents(pbmc) <- 'bulk'
p4 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(30) 
Idents(pbmc) <- 'bulk.all'
p5 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
p2 + p3 + p4 +p5 

Idents(pbmc) <- 'singler.fine'
DotPlot(pbmc , features = c("PDGFRA", "LUM", "DCN", "VIM",  "COL1A2", "PMEL", "MLANA", "TYRP1",  "DCT"))
Idents(pbmc) <- 'bulk'
DotPlot(pbmc , features = c("PDGFRA", "LUM", "DCN", "VIM",  "COL1A2", "PMEL", "MLANA", "TYRP1",  "DCT"))
Idents(pbmc) <- 'bulk.all'
DotPlot(pbmc , features = c("PDGFRA", "LUM", "DCN", "VIM",  "COL1A2", "PMEL", "MLANA", "TYRP1",  "DCT"))


Idents(pbmc) <- 'RNA_snn_res.0.3'
clusterX <- pbmc[,grep("^3|^5", Idents(pbmc))]
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(15) 
Idents(clusterX) <- 'RNA_snn_res.0.3'
p2 <- DimPlot(clusterX, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(38) 
Idents(clusterX) <- 'singler.fine'
p3 <- DimPlot(clusterX, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(29) 
Idents(clusterX) <- 'bulk'
p4 <- DimPlot(clusterX, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(30) 
Idents(clusterX) <- 'bulk.all'
p5 <- DimPlot(clusterX, cols = pal.celltype, label = T) + NoLegend()
p2 + p3 + p4 +p5 

table(clusterX$singler.fine, clusterX$bulk)
pred.df <- data.frame(SingleR = clusterX$singler.fine, bulk = clusterX$bulk, 
                      freq = c(1))
#sankey plot
library(ggalluvial)
sankeyp <- ggplot(data = pred.df,
                  aes(axis1 = SingleR, axis2 = bulk, y = freq)) +
  geom_alluvium(aes(fill = bulk)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("SingleR", "bulk"),
                   expand = c(0.15, 0.05)) +
  theme_void()







#pseudobulk 
Idents(clusterX) <- 'orig.ident'
pseudobulk <- AggregateExpression(clusterX, assays = "RNA")
pseudobulk <- round(scale_to_million(pseudobulk$RNA),4)
decon.pseudobulk <- cibersort(sig.mat[rownames(meanGini)[1:2000],grep("Plasma",colnames(sig.mat), invert = T)], as.matrix(pseudobulk))
decon.pseudobulk <- t(decon.pseudobulk[,1:(ncol(decon.pseudobulk)-3)])
decon.pseudobulk



####
ref.bulk <- log2(sig.mat[rownames(meanGini)[1:1000],]+1)
keep <- intersect(rownames(pbmc.data), rownames(ref.bulk))
pbmc.decon <- pbmc.data[keep, ]
decon.cibersort <- cibersort(ref.bulk, as.matrix(pbmc.decon))
####








###ENCODE data only########
rownames(target_tissue) <- gsub("-", "_", rownames(target_tissue))
target_tissue <- target_tissue[grep("^MT_", rownames(target_tissue), invert = T),]
target_tissue <- target_tissue[grep("^RPL[0-9]", rownames(target_tissue), invert = T),]
target_tissue <- target_tissue[grep("^MIR[0-9]", rownames(target_tissue), invert = T),]
target_tissue <- target_tissue[grep("^LINC[0-9]", rownames(target_tissue), invert = T),]
target_tissue <- target_tissue[grep("^RNU[0-9]*", rownames(target_tissue), invert = T),]
target_tissue <- target_tissue[grep("\\.", rownames(target_tissue), invert = T),]
topVarGenes <- head(order(rowVars(as.matrix(target_tissue)), decreasing = TRUE), 5000)
rownames(target_tissue)[topVarGenes]
data <- as.data.frame(t(target_tissue[topVarGenes,]))
data$CellType <- pdataENCODE$cell
data$CellType <- as.factor(data$CellType)
table(data$CellType)

set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T) 
print(rf)
rf
varImpPlot(rf,
           sort = T,
           n.var = 30,
           main = "Top 30 - Variable Importance")
meanGini <- importance(rf)
meanGini <- meanGini[order(meanGini[,"MeanDecreaseAccuracy"], decreasing = T),]
meanGini2 <- importance(rf)[order(importance(rf)[,"MeanDecreaseGini"], decreasing = T),]

####sig.mat ######
sig.mat <- target_tissue[, ]
pdataENCODE <- cbind(pdataENCODE, t(sig.mat))
sig.mat <- aggregate(x = pdataENCODE[,c(5:ncol(pdataENCODE))],by = list(pdataENCODE$cell), FUN = median)
rownames(sig.mat) <- sig.mat$Group.1
sig.mat <- sig.mat[,-1]
sig.mat <- t(sig.mat)


###in silico #####
library(CIBERSORT)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(31) 

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:500],]), rownames(target_tissue))
tpm.decon <- target_tissue[keep,]
cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:500],], as.matrix(tpm.decon))
cibersort.df <- t(round(cibersort.df[,1:(ncol(cibersort.df)-3)], 4))

library(pheatmap)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(31)  
colnames(cibersort.df) <- pdataENCODE$id
mat_col <- data.frame(group = pdataENCODE$cell)
rownames(mat_col) <- colnames(cibersort.df)
mat_colors <- list(group = pal.celltype)
names(mat_colors$group) <- unique(pdataENCODE$cell)
#cibersort.df <- round(cibersort.df*100)
pheatmap(cibersort.df[order(rownames(cibersort.df)),order(pdataENCODE$cell)], 
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = "black", 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "RF500 CIBERSORT ENCODE" )


