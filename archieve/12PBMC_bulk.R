library(EnsDb.Hsapiens.v86)
library(pheatmap)

pbmc <- read.table("./dataset/COVID/GSE107011_Processed_data_TPM.txt.gz", sep = "\t", row.names = 1, header = T)
pbmc$ensg <- substr(rownames(pbmc), start = 1L, stop = 15L)
pbmc <- aggregate(x = pbmc[,-ncol(pbmc)],by = list(pbmc$ensg), FUN = median)
rownames(pbmc) <- pbmc$Group.1
pbmc <- pbmc[,-1]


ensembl.genes <-  rownames(pbmc)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(pbmc))
pbmc <- pbmc[keep,]
data.frame(rownames(pbmc), geneIDs$GENEID)
pbmc$symbol <- geneIDs$SYMBOL
pbmc <- aggregate(x = pbmc[,-ncol(pbmc)],by = list(pbmc$symbol), FUN = median)
rownames(pbmc) <- pbmc$Group.1
pbmc <- pbmc[,-1]
pbmc[1:10,1:10]

goi <- c("DAPP1", "CST3", "FGL2", "GCH1", "CIITA", "UPP1",  "RN7SL1")
pheatmap(pbmc[goi,], scale = "row", 
         cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","white","darkred"))(100), 
         border_color      = NA,
         show_colnames     = T,show_rownames     = T,
         main = "PBMC" )


library(GEOquery)
GSEset <- getGEO("GSE107011", destdir = "./dataset/COVID/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 40:42)]
sampleinfo$title
colnames(sampleinfo)[3] <- "celltype"
table(sampleinfo$celltype) #n = 30 cell types

save(pbmc,sampleinfo, file = "./dataset/COVID/GSE107011.RData")

# Data frame with column annotations.
load(file = "./dataset/COVID/GSE107011.RData")

mat_col <- data.frame(group = sampleinfo$celltype)
rownames(mat_col) <- colnames(pbmc)
# List with colors for each annotation.
mat_colors <- list(group = colorRampPalette(brewer.pal(8, "Set3"))(30) )
names(mat_colors$group) <- unique(sampleinfo$celltype)
pheatmap(pbmc[goi,], scale = "row", 
         cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","white","darkred"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )
t(pbmc[goi,])
sampleinfo <- cbind(sampleinfo, t(pbmc[goi,]))

tmp <- sampleinfo[grep("PBMC|B cells|monocytes|dendritic", sampleinfo$celltype),]
boxplot(tmp$DAPP1 ~ tmp$celltype)






###1. Random forest , marker selection####
load(file = "./dataset/COVID/PBMC_GEO.RData")
library(randomForest)
library(caret)
library(edgeR)
library(matrixStats)
library(EnsDb.Hsapiens.v86)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "SYMBOL", columns = c("SYMBOL", "GENEBIOTYPE"))
geneIDs <- geneIDs[grep("protein_coding", geneIDs$GENEBIOTYPE),]
keep <- intersect(geneIDs$SYMBOL, rownames(tpm))
rownames(geneIDs) <- geneIDs$SYMBOL
geneIDs <- geneIDs[keep, ]
tpm <- tpm[keep,]

tpm <- tpm[grep("\\.", rownames(tpm), invert = T),]
tpm <- tpm[grep("^MT-", rownames(tpm), invert = T),]
tpm <- tpm[grep("^HLA-", rownames(tpm), invert = T),]
tpm <- tpm[grep("-AS", rownames(tpm), invert = T),]
tpm <- tpm[grep("PALM2-AKAP2", rownames(tpm), invert = T),]
tpm <- tpm[grep("KLRC4-KLRK1", rownames(tpm), invert = T),]
topVarGenes <- head(order(rowVars(as.matrix(tpm)), decreasing = TRUE), 5000)
rownames(tpm)[topVarGenes]
data <- as.data.frame(t(tpm[topVarGenes,]))
data$CellType <- pdata$cell
data <- data[grep("PBMC|Th1Th17", data$CellType, invert = T),]
data$CellType <- as.factor(data$CellType)
table(data$CellType)

#set.seed(222)
#ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
#train <- data[ind==1,]
#test <- data[ind==2,]

#http://www.ehbio.com/ML/randomForest.html
set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T,  mtry = 157) 
print(rf)
rf
set.seed(304)
tuneRF(data[,-ncol(data)], data$CellType, ntreeTry=500, stepFactor=1.5, improve=1e-5)
set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T,  mtry = 157) 
print(rf)

#p1 <- predict(rf, train)
#confusionMatrix(p1, train$ Species)
plot(rf)
t <- tuneRF(train[,-5], train[,5],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)
hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")

varImpPlot(rf,
           sort = T,
           n.var = 100,
           main = "Top 20 - Variable Importance")
meanGini <- importance(rf)
meanGini <- meanGini[order(meanGini[,"MeanDecreaseGini"], decreasing = T),]
MDSplot(rf, data$CellType)
RF.genelist1000 <- rownames(meanGini)[1:1000]
RF.genelist500 <- rownames(meanGini)[1:500]
RF.genelist100 <- rownames(meanGini)[1:100]
RF.genelist157 <- rownames(meanGini)[1:157]

save(RF.genelist100, RF.genelist1000, RF.genelist157, RF.genelist500, file = "./dataset/COVID/RF_sig.RData")

###2. cut off ####
#load(file = "./PBMC_GEO.RData")
library(EnsDb.Hsapiens.v86)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "SYMBOL", columns = c("SYMBOL", "GENEBIOTYPE"))
geneIDs <- geneIDs[grep("protein_coding", geneIDs$GENEBIOTYPE),]
keep <- intersect(geneIDs$SYMBOL, rownames(tpm))
rownames(geneIDs) <- geneIDs$SYMBOL
geneIDs <- geneIDs[keep, ]
tpm <- tpm[keep,]

idx <- grep("PBMC|Th1Th17", pdata$cell, invert = T)
library(RColorBrewer)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(27)  
  
plotMDS(log2(tpm[,idx]+1), top=2000, gene.selection="common", pch = 19,
        col=pal.celltype[factor(pdata[idx,]$cell)])
legend("top", legend=levels(factor(pdata[idx,]$cell)), text.col=pal.celltype,
       bg="white", cex=0.7)

tpm <- tpm[,idx]
pdata <- pdata[idx,]

cutoff.genelist <- c()
for (i in unique(pdata$cell)) {
  print(i)
  idx <- grep(i, pdata$cell)
  idx2 <- grep(i, pdata$cell, invert = T)
  #tmp.list <- rownames(tpm)[which(rowMeans( tpm[,idx])- rowMeans( tpm[,idx2]) > 100 &  rowMeans( tpm[,idx2]) < 50) ]
  tmp.list <- rownames(tpm)[order(rowMeans( tpm[,idx])/ rowMeans( tpm[,idx2]) , decreasing = T)]
  cutoff.genelist <- c(cutoff.genelist, tmp.list[1:100])
  #cutoff.genelist <- c(cutoff.genelist, tmp.list)
}
cutoff.genelist <- unique(cutoff.genelist)

genelist <- unique(c(cutoff.genelist, RF.genelist))
intersect(genelist, rownames(abis_sig))

save(cutoff.genelist, RF.genelist, genelist, file = "./dataset/COVID/PBMC_GEO_signature.RData")


###3. marker matrix ####
plotMDS(log2(tpm[genelist,]+1), top=1000, gene.selection="common", pch = 19,
        col=pal.celltype[factor(pdata$cell)])
plotMDS(log2(tpm[RF.genelist,]+1), top=1000, gene.selection="common", pch = 19,
        col=pal.celltype[factor(pdata$cell)])
plotMDS(log2(tpm[cutoff.genelist,]+1), top=1000, gene.selection="common", pch = 19,
        col=pal.celltype[factor(pdata$cell)])

sig.mat <- data.frame(mock = rep(1, nrow(tpm)))
for (i in unique(pdata$cell)) {
  print(i)
  idx <- grep(i, pdata$cell)
  sig.mat <- cbind(sig.mat,  rowMeans( tpm[,idx]) )
}
#sig.mat <- sig.mat[genelist,]
sig.mat <- sig.mat[,-1]
colnames(sig.mat) <- unique(pdata$cell)

save(sig.mat, file = "./dataset/COVID/PBMC_GEO_sigMat.rds")


###4. refine the markers (add Th1)####
sig.mat <- sig.mat[RF.genelist157,]
sig.mat <- sig.mat[,grep("PBMC|Th1Th17", colnames(sig.mat), invert = T)]

cell.name <- c()
for (i in 1:nrow(sig.mat)) {
  cell.name <- c(cell.name,names(which.max(sig.mat[i,])))
}
cell.manifest <- data.frame(gene = rownames(sig.mat), cellname = cell.name)
table(cell.manifest$cellname)

