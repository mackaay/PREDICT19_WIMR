
library(randomForest)
library(caret)
library(e1071)
library(sigFeature)
library(matrixStats)
library(CIBERSORT)
library(RColorBrewer)

load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/all.data.RData")
load(file = "./dataset/COVID/RF_meanGini.rds")
load(file = "./dataset/COVID/RF_meanGini2.rds")
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))

load(file = "./dataset/COVID/ABIS_RF_sigMat_adjPMN.rds")
load(file = "./dataset/COVID/all.adjusted_PMN.RData")
all.tpm <- adjusted.tpm
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
load(file = "./dataset/COVID/RF_meanGini_adjPMN.rds")
load(file = "./dataset/COVID/RF_meanGini2_adjPMN.rds")




load(file = "./dataset/COVID/all_adjusted.RData")
load(file = "./dataset/COVID/RF_meanGini_adj.rds")
load( file = "./dataset/COVID/SigMatrix_adj.rds")
all.tpm <- adjusted.tpm
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))

load(file = "./dataset/COVID/pal.celltype.rds")




rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
all.tpm <- all.tpm[grep("^MT_", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^RPL[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^MIR[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^LINC[0-9]", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("^RNU[0-9]*", rownames(all.tpm), invert = T),]
all.tpm <- all.tpm[grep("\\.", rownames(all.tpm), invert = T),]
#topVarGenes <- head(order(rowVars(as.matrix(all.tpm)), decreasing = TRUE), 5000)
#rownames(all.tpm)[topVarGenes]
#data <- all.tpm[topVarGenes,]
#data <- as.data.frame(t(data))
#data$CellType <- all.pdata$cell
#data$CellType <- as.factor(data$CellType)


###iteration #####
parallelly::supportsMulticore()
library(future)
plan(multicore)
accuracy <- data.frame(ngene =c(), Cor = c(), RMSE = c())
for (i in 2:250) {
  keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:(i*10)],]), rownames(all.tpm))
  tpm.decon <- all.tpm[keep,]
  cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:(i*10)],], as.matrix(tpm.decon), QN = F)
  print((i*10))
  print(median(cibersort.df[,"Correlation"]))
  print(median(cibersort.df[,"RMSE"]) )
  tmp.acc <- data.frame(ngene = rep((i*10), nrow(cibersort.df)), 
                        Cor = cibersort.df[,"Correlation"], 
                        RMSE = cibersort.df[,"RMSE"])
  accuracy <- rbind(accuracy, tmp.acc)
  #save(accuracy, file = "./dataset/COVID/SVM_adjPMN_accuracy.rds")
  save(accuracy, file = "./dataset/COVID/SVM_accuracy_adj_noQN.rds")
}

tail(accuracy)

for (i in 197:250) {
  keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:(i*10)],]), rownames(all.tpm))
  tpm.decon <- all.tpm[keep,]
  cibersort.df  <- cibersort(sig.mat[rownames(meanGini)[1:(i*10)],], as.matrix(tpm.decon), QN = F)
  print((i*10))
  print(median(cibersort.df[,"Correlation"]))
  print(median(cibersort.df[,"RMSE"]) )
  tmp.acc <- data.frame(ngene = rep((i*10), nrow(cibersort.df)), 
                        Cor = cibersort.df[,"Correlation"], 
                        RMSE = cibersort.df[,"RMSE"])
  accuracy <- rbind(accuracy, tmp.acc)
  #save(accuracy, file = "./dataset/COVID/SVM_adjPMN_accuracy.rds")
  save(accuracy, file = "./dataset/COVID/SVM_accuracy_adj_noQN_3.rds")
}

#merge three SVM accuracy files together
load(file = "./dataset/COVID/SVM_accuracy_adj_noQN.rds")
accuracy1 <- accuracy
load(file = "./dataset/COVID/SVM_accuracy_adj_noQN_2.rds")
accuracy2 <- accuracy
load(file = "./dataset/COVID/SVM_accuracy_adj_noQN_3.rds")
accuracy3 <- accuracy

accuracy <- rbind(accuracy1, accuracy2)
accuracy <- rbind(accuracy, accuracy3)

save(accuracy, file = "./dataset/COVID/SVM_accuracy_adj_noQN_20250120.rds")
write.csv(accuracy, file = "./dataset/COVID/output/12accuracy.csv")


###20241113plot####
load(file = "./dataset/COVID/SVM_accuracy.rds") #20:2000
load("./dataset/COVID/SVM_adjPMN_accuracy.rds") #20:2610
load("./dataset/COVID/SVM_accuracy_noQN.rds") #20:2500
load(file = "./dataset/COVID/SVM_accuracy_adj_noQN.rds")
load(file = "./dataset/COVID/SVM_accuracy_adj_noQN_20250120.rds")


library(ggplot2)
library(dplyr)
accuracy$CellType <- rep(all.pdata$cell, nrow(accuracy)/nrow(all.pdata))
accuracy[nrow(accuracy),]
accuracy.df <- accuracy %>% 
  group_by(ngene) %>% 
  summarise(Correlation = median(Cor, na.rm = TRUE), RMSE= median(RMSE, na.rm = TRUE))
ggplot(data = accuracy.df, aes(x = ngene, y = Correlation))+
  geom_line()+ ylim(0.2,1)+ 
  geom_vline(xintercept=c(2280), linetype="dotted", color = "darkred")+
  labs(x = "No. of Genes")
  #geom_vline(xintercept=c(325),  color = "yellow", size = 25, alpha = 0.2)
ggplot(data = accuracy.df, aes(x = ngene, y = RMSE))+
  geom_line()+ ylim(0.25,1)+ 
  geom_vline(xintercept=c(2280), linetype="dotted", color = "darkred")+
  labs(x = "No. of Genes")
dev.off()


accuracy.df <- accuracy %>%
  group_by(ngene, CellType) %>%
  summarise(Correlation = median(Cor, na.rm = TRUE), RMSE = median(RMSE, na.rm = TRUE))
accuracy.cutoff <- accuracy.df %>% 
  group_by(ngene) %>%
  summarise(Correlation = median(Correlation, na.rm = TRUE), RMSE = median(RMSE, na.rm = TRUE))
accuracy.cutoff[order(accuracy.cutoff$Correlation, decreasing = T),]
accuracy.cutoff[order(accuracy.cutoff$RMSE, decreasing = F),]

ggplot(data = accuracy.df, aes(x = ngene, y = Correlation))+
  geom_line()+ ylim(0.5,1) + 
  geom_vline(xintercept=c(2280), linetype="dotted", color = "darkred")+
  facet_wrap(.~CellType,  ncol = 7 )+
  labs(x = "No. of Genes")
ggplot(data = accuracy.df, aes(x = ngene, y = RMSE))+
  geom_line()+ ylim(0,1) + 
  geom_vline(xintercept=c(2280), linetype="dotted", color = "darkred")+
  facet_wrap(.~CellType )+
  labs(x = "No. of Genes")
ggplot(data = accuracy.df, aes(x = ngene, y = Correlation, group = CellType, color = CellType))+
  geom_line(aes(color = CellType))+ ylim(0.2,1) + 
  geom_vline(xintercept=c(2280), linetype="dotted", color = "darkred")+
  scale_color_manual(values = pal.celltype)
ggplot(data = accuracy.df, aes(x = ngene, y = RMSE, group = CellType, color = CellType))+
  geom_line(aes(color = CellType))+ ylim(0,1) + 
  geom_vline(xintercept=c(2280), linetype="dotted", color = "darkred")+
  scale_color_manual(values = pal.celltype)







###intersection of each pairwise comparision#####
unique(all.pdata$cell)
idx <- grep("CD8Tnaive", all.pdata$cell)
x <- t(data[,])
y <- rep(-1, nrow(x))
y[idx] = 1
pvals <- sigFeaturePvalue(x,y)
hist(unlist(pvals),col="skyblue", #breaks=seq(0,0.08,0.0015),
     xlab="p value",ylab="Frequency",main="")
temp.sig <- rownames(data)[which(unlist(pvals) < 0.05)]
sig.feature <- temp.sig 
print(length(sig.feature))
for (i in unique(all.pdata$cell)[2:length(unique(all.pdata$cell))]) {
  idx <- grep(i, all.pdata$cell)
  x <- t(data[,])
  y <- rep(-1, nrow(x))
  y[idx] = 1
  pvals <- sigFeaturePvalue(x,y)
  hist(unlist(pvals),col="skyblue", #breaks=seq(0,0.08,0.0015),
       xlab="p value",ylab="Frequency",main="")
  temp.sig <- rownames(data)[which(unlist(pvals) < 0.05)]
  sig.feature <- intersect(sig.feature, temp.sig)
  print(length(sig.feature))
}
sig.feature


###merge of each pairwise comparision#####
unique(all.pdata$cell)
idx <- grep("CD8Tnaive", all.pdata$cell)
x <- t(data[,])
y <- rep(-1, nrow(x))
y[idx] = 1
pvals <- sigFeaturePvalue(x,y)
hist(unlist(pvals),col="skyblue", #breaks=seq(0,0.08,0.0015),
     xlab="p value",ylab="Frequency",main="")
temp.sig <- rownames(data)[which(unlist(pvals) < 0.0001)]
sig.feature <- temp.sig 
print(length(sig.feature))
for (i in unique(all.pdata$cell)[2:length(unique(all.pdata$cell))]) {
  idx <- grep(i, all.pdata$cell)
  x <- t(data[,])
  y <- rep(-1, nrow(x))
  y[idx] = 1
  pvals <- sigFeaturePvalue(x,y)
  hist(unlist(pvals),col="skyblue", #breaks=seq(0,0.08,0.0015),
       xlab="p value",ylab="Frequency",main="")
  temp.sig <- rownames(data)[which(unlist(pvals) < 0.001)]
  sig.feature <- unique(c(sig.feature, temp.sig))
  print(length(sig.feature))
}
sig.feature


###Embedded methods LASSO ####
library(glmnet)
data <- all.tpm[,]
data <- log2(data +1)
#x <- scale(data )
x <- t(x)
set.seed(2)
lasso_mod = glmnet(x = x, 
                   y = all.pdata$cell, 
                   #lambda = grid,
                   family = "multinomial",
                   alpha = 1 
) # Fit lasso model on training data
plot(lasso_mod, xvar = "lambda", label = T)    # Draw plot of coefficients
print(lasso_mod)
coeff <- coef(lasso_mod)
which(abs(coeff[[1]][,10]) > 0)
names(coeff[[1]][,10])[-which(coeff[[1]][,10] == 0)]
varImp(lasso_mod)

library(glmnet)
data <- all.tpm[,]
data <- log2(data +1)
x <- scale(data )
x <- t(x)

sig.feature <- c()
for (i in unique(all.pdata$cell)[1:length(unique(all.pdata$cell))]) {
  idx <- grep(i, all.pdata$cell)
  y <- rep(-1, nrow(x))
  y[idx] = 1
  set.seed(2)
  lasso_mod = glmnet(x = x, 
                     y = y, 
                     #lambda = grid,
                     family = "binomial",
                     alpha = 1 
  ) # Fit lasso model on training data
  #plot(lasso_mod, xvar = "lambda", label = T)    # Draw plot of coefficients
  #print(lasso_mod)
  sig.feature <- c(names(coef(lasso_mod)[,10])[which(abs(coef(lasso_mod)[,10]) > 0)], sig.feature)
  print(length(sig.feature))
}



save(sig.feature, file = "./dataset/COVID/SVM_sigfeature.rds")



system.time(sigFeature(x[,which(unlist(pvals) < 0.0001)], y))
sigfeatureRankedList <- sigFeature(x[,which(unlist(pvals) < 0.0001)], y)
