###1. Random forest , marker selection####
load(file = "./dataset/COVID/all_adjusted.RData")
load(file = "./dataset/COVID/pal.celltype.rds")
all.tpm <- adjusted.tpm
write.csv(all.pdata, file = "./plots/all.pdata.csv")




library(randomForest)
library(caret)
library(edgeR)
library(matrixStats)

rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
all.tpm <- all.tpm[grep("^MT_", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^RPL[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^MIR[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^LINC[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^RNU[0-9]*", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("\\.", rownames(all.tpm), invert = T),]
#topVarGenes <- head(order(rowVars(as.matrix(all.tpm)), decreasing = TRUE), 10000)
#rownames(all.tpm)[topVarGenes]
data <- as.data.frame(t(all.tpm[1:10000,]))
data$CellType <- all.pdata$cell
data$CellType <- as.factor(data$CellType)
table(data$CellType)



set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T) 
print(rf)
rf
meanGini <- importance(rf)
meanGini <- meanGini[order(meanGini[,"MeanDecreaseAccuracy"], decreasing = T),]
meangini.name1 <-rownames(meanGini)[1:5000]


data <- as.data.frame(t(all.tpm[10001:nrow(all.tpm),]))
data$CellType <- all.pdata$cell
data$CellType <- as.factor(data$CellType)
table(data$CellType)
set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T) 
meanGini <- importance(rf)
meanGini <- meanGini[order(meanGini[,"MeanDecreaseAccuracy"], decreasing = T),]
meangini.name2 <-rownames(meanGini)[1:((nrow(all.tpm)-10000)/2)]

#merge 1 and 2
data <- as.data.frame(t(all.tpm[c(meangini.name1, meangini.name2),]))
data$CellType <- all.pdata$cell
data$CellType <- as.factor(data$CellType)
table(data$CellType)
set.seed(304)
rf <- randomForest(CellType~., data=data, proximity=TRUE, importance = T) 
meanGini <- importance(rf)
meanGini <- meanGini[order(meanGini[,"MeanDecreaseAccuracy"], decreasing = T),]
meanGini2 <- importance(rf)[order(importance(rf)[,"MeanDecreaseGini"], decreasing = T),]



save(meanGini, file = "./dataset/COVID/RF_meanGini.rds")
save(meanGini2, file = "./dataset/COVID/RF_meanGini2.rds")
save(meanGini, file = "./dataset/COVID/RF_meanGini_adjPMN.rds")
save(meanGini2, file = "./dataset/COVID/RF_meanGini2_adjPMN.rds")
save(meanGini, file = "./dataset/COVID/RF_meanGini_adj.rds")
save(meanGini2, file = "./dataset/COVID/RF_meanGini2_adj.rds")







###Signature Matrix #####
sig.mat <- all.tpm[, ]
all.pdata <- cbind(all.pdata, t(sig.mat))
sig.mat <- aggregate(x = all.pdata[,c(8:ncol(all.pdata))],by = list(all.pdata$cell), FUN = median)
rownames(sig.mat) <- sig.mat$Group.1
sig.mat <- sig.mat[,-1]
sig.mat <- t(sig.mat)

save(sig.mat, file = "./dataset/COVID/ABIS_RF_sigMat.rds")
save(sig.mat, file = "./dataset/COVID/ABIS_RF_sigMat_adjPMN.rds")
save(sig.mat, file = "./dataset/COVID/SigMatrix_adj.rds")

topVarGenes <- head(order(rowVars(as.matrix(sig.mat)), decreasing = TRUE), 10000)
topVarGenes <- rownames(sig.mat)[topVarGenes]
save(topVarGenes, file = "./dataset/COVID/topVarGenes.rds")



