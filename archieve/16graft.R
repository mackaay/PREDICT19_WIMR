library(pheatmap)
library(RColorBrewer)
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/gene_length.rds")
load(file = "./dataset/COVID/RF_meanGini.rds")
load(file = "./dataset/COVID/RF_meanGini2.rds")
niteration <- 1220




load(file = "./dataset/COVID/all_adjusted.RData")
load(file = "./dataset/COVID/RF_meanGini_adj.rds")
load( file = "./dataset/COVID/SigMatrix_adj.rds")
all.tpm <- adjusted.tpm
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))



load(file = "./dataset/COVID/all_adjusted.RData")
load(file = "./dataset/COVID/RF_meanGini_adj.rds")
load( file = "./dataset/COVID/SigMatrix_adj.rds")
all.tpm <- adjusted.tpm
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
niteration <- 2280
load(file = "./dataset/COVID/pal.celltype.rds")



##KidneyGraft_GSE248752####
library(CIBERSORT)
load(file = "./dataset/COVID/KidneyGraft_GSE248752/KidneyGraft_GSE248752.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

#PBMC or whole blood??
#sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

mat_col <- data.frame(Treatment = pdata$`treatment:ch1`)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Treatment = brewer.pal(3, "Set1")
                   )
names(mat_colors$Treatment) <- unique(pdata$`treatment:ch1`)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.df[, order(pdata$treatment)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = " KidneyGraft_GSE248752" )

pdata <- cbind(pdata, t(decon.df))

pdata.df <- pdata[,3:ncol(pdata)]
pdata.df <- reshape2::melt(pdata.df)
colnames(pdata.df) <- c("treatment", "CellType", "Proportion")

library(ggplot2)
library(dplyr)
library(ggpubr)
pdata.df[complete.cases(pdata.df),]%>%
  ggplot( aes(x=treatment, y=Proportion, fill=treatment)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("CellType")+
  geom_signif(comparisons = list(c("Nill", "Rej")), y_position = c(0.35))



##KidneyGraft_GSE120396####
library(CIBERSORT)
load(file = "./dataset/COVID/KidneyGraft_GSE120396/KidneyGraft_GSE120396.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

#PBMC or whole blood??
#sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

mat_col <- data.frame(Treatment = pdata$characteristics_ch1.1)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Treatment = brewer.pal(3, "Set1")
)
names(mat_colors$Treatment) <- unique(pdata$characteristics_ch1.1)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.df[, order(pdata$characteristics_ch1.1)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "KidneyGraft_GSE120396" )

pdata <- cbind(pdata, t(decon.df))

pdata.df <- pdata[,3:ncol(pdata)]
pdata.df <- reshape2::melt(pdata.df)
colnames(pdata.df) <- c("ACR_3m", "CellType", "Proportion")
pdata.df$ACR_3m <- gsub("acr at 3m: ", "", pdata.df$ACR_3m)

library(ggplot2)
library(dplyr)
library(ggpubr)
pdata.df[complete.cases(pdata.df),]%>%
  ggplot( aes(x=ACR_3m, y=Proportion, fill=ACR_3m)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("CellType")+
  geom_signif(comparisons = list(c("No", "Yes")), y_position = c(0.3))



##KidneyGraft_GSE120649####
library(CIBERSORT)
load(file = "./dataset/COVID/KidneyGraft_GSE120649/KidneyGraft_GSE120649.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))
pdata$group <- c(rep("AntibodyRej", 6), rep("Stable", 6), rep("TcellRej",4))

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

#PBMC or whole blood??
#sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

mat_col <- data.frame(Treatment = pdata$group)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Treatment = brewer.pal(3, "Set1")
)
names(mat_colors$Treatment) <- unique(pdata$group)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.df[, order(pdata$group)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "KidneyGraft_GSE120396" )

pdata <- cbind(pdata, t(decon.df))

pdata.df <- pdata[,3:ncol(pdata)]
pdata.df <- reshape2::melt(pdata.df)
colnames(pdata.df) <- c("group", "CellType", "Proportion")

library(ggplot2)
library(dplyr)
library(ggpubr)
pdata.df[complete.cases(pdata.df),]%>%
  ggplot( aes(x=group, y=Proportion, fill=group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  scale_fill_manual(values=c("#85C1E9", "#F7DC6F", "firebrick"))+
  facet_wrap("CellType")+
  geom_signif(comparisons = list(c("AntibodyRej", "Stable")), y_position = c(0.21))





##KidneyGraft_GSE175718####
library(CIBERSORT)
load( file = "./dataset/COVID/KidneyGraft_GSE175718/KidneyGraft_GSE175718.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))
pdata$`abmrh:ch1` <- as.numeric(pdata$`abmrh:ch1`)
pdata$`tcmr:ch1` <- as.numeric(pdata$`tcmr:ch1`)
pdata$Rej <- "No"
pdata$Rej[pdata$`abmrh:ch1` == 1] <- "AbMR"
pdata$Rej[grep(1, pdata$`tcmr:ch1`)] <- "TCMR"


library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

#PBMC or whole blood??
#sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

mat_col <- data.frame(Treatment = pdata$Rej)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Treatment = brewer.pal(3, "Set1")
)
names(mat_colors$Treatment) <- unique(pdata$Rej)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.df[, order(pdata$Rej)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = " KidneyGraft_GSE248752" )

pdata <- cbind(pdata, t(decon.df))

pdata.df <- pdata[,8:ncol(pdata)]
pdata.df <- reshape2::melt(pdata.df)
colnames(pdata.df) <- c("treatment", "CellType", "Proportion")

library(ggplot2)
library(dplyr)
library(ggpubr)
pdata.df[complete.cases(pdata.df),]%>%
  ggplot( aes(x=treatment, y=Proportion, fill=treatment)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  scale_fill_manual(values=c("#85C1E9", "#F7DC6F", "firebrick"))+
  facet_wrap("CellType")+
  geom_signif(comparisons = list(c("AbMR", "No")), y_position = c(0.35))








##PBMC night shift Nutrition GSE225493####
library(CIBERSORT)
load(file = "./dataset/COVID/PBMC_NightShiftNutrition_GSE225493/PBMC_NightShiftNutrition_GSE225493.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )


sig.mat <- sig.mat[,grep("Neutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

mat_col <- data.frame(Treatment = pdata$`treatment:ch1`,
                      Timepoint = pdata$`timepoint:ch1`)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Treatment = brewer.pal(3, "Set1"),
                   Timepoint = brewer.pal(3, "Dark2"))
names(mat_colors$Treatment) <- unique(pdata$`treatment:ch1`)
names(mat_colors$Timepoint) <- unique(pdata$`timepoint:ch1`)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

colnames(pdata)[4:5] <- c("timepoint", "treatment")
pheatmap(decon.1.df[, order(pdata$treatment)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "PBMC night shift Nutrition CIBERSORT meanGini 1220" )

library(ggplot2)
library(ggpubr)
pdata <- cbind(pdata, decon.cibersort)
colnames(pdata)[4:5] <- c("timepoint", "treatment")

pdata$mono <- pdata$Mono_C+ pdata$Mono_I + pdata$Mono_NC
pdata[complete.cases(pdata),]%>%
  ggplot( aes(x=treatment, y=Mono_C, fill=treatment)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_C ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F", "firebrick"))+
  facet_wrap("timepoint")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Control", "Intervention")), y_position = c(0.35))
pdata[complete.cases(pdata),]%>%
  ggplot( aes(x=treatment, y=CD4Tnaive, fill=treatment)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("CD4Tnaive ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F", "firebrick"))+
  facet_wrap("timepoint")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Control", "Intervention")), y_position = c(0.20))
pdata[complete.cases(pdata),]%>%
  ggplot( aes(x=treatment, y=Th2, fill=treatment)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Th2 ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F", "firebrick"))+
  facet_wrap("timepoint")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Control", "Intervention")), y_position = c(0.25))



library(reshape2)
pdata <- melt(pdata[,3:ncol(pdata)])
colnames(pdata)[1] <- "SID" 
for (i in 1:length(unique(pdata$variable))-3) {
  p <-   ggplot(pdata[grep(paste0(unique(pdata$variable)[i], "$"), pdata$variable),], aes(x=timepoint, y=value, group=treatment, colour=treatment ) ) + 
    geom_line(size=1)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F"))+ 
    geom_point(size=1.5)+
    facet_wrap(SID~.) +
    theme(axis.text.x = element_text(angle = 90)) + ylab(unique(pdata$variable)[i])
  print(p)
}


ggplot(pdata[grep("Th1$", pdata$variable),], aes(x=timepoint, y=value, group=treatment, colour=treatment ) ) + 
  geom_line(size=1)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F"))+ 
  geom_point(size=1.5)+
  facet_wrap(SID~.) +
  theme(axis.text.x = element_text(angle = 90)) 

ggplot(pdata[grep("Mono_C", pdata$variable),], aes(x=timepoint, y=value, group=treatment, colour=treatment ) ) + 
  geom_line(size=1)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F"))+ 
  geom_point(size=1.5)+
  facet_wrap(SID~.) +
  theme(axis.text.x = element_text(angle = 90)) 

ggplot(pdata[grep("CD8Tcm", pdata$variable),], aes(x=timepoint, y=value, group=treatment, colour=treatment ) ) + 
  geom_line(size=1)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F"))+ 
  geom_point(size=1.5)+
  facet_wrap(SID~.) +
  theme(axis.text.x = element_text(angle = 90)) 

ggplot(pdata[grep("CD8Teff", pdata$variable),], aes(x=timepoint, y=value, group=treatment, colour=treatment ) ) + 
  geom_line(size=1)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F"))+ 
  geom_point(size=1.5)+
  facet_wrap(SID~.) +
  theme(axis.text.x = element_text(angle = 90)) 





##new-onset T1D GSE124400####
library(CIBERSORT)
load(file = "./dataset/COVID/T1D_GSE124400/T1D_GSE124400.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))
pdata

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

mat_col <- data.frame(Race = pdata$`race:ch1`
                      )
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Race = colorRampPalette( brewer.pal(8, "Set1"))(10)
                   )
names(mat_colors$Race) <- unique(pdata$`race:ch1`)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.df[, order(pdata$`rate of c-peptide change:ch1`)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "T1D_GSE124400 CIBERSORT meanGini 1220" )





##T1D_MS_GSE60424####
library(CIBERSORT)
load(file = "./dataset/COVID/T1D_MS_GSE60424/T1D_MS_GSE60424.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))
pdata

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

mat_col <- data.frame(CellType = pdata$`celltype:ch1`, 
                      Disease = pdata$`diseasestatus:ch1`
)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(CellType = colorRampPalette( brewer.pal(8, "Accent"))(7), 
                   Disease = colorRampPalette( brewer.pal(8, "Set1"))(6)
)
names(mat_colors$CellType) <- unique( pdata$`celltype:ch1`)
names(mat_colors$Disease) <- unique( pdata$`diseasestatus:ch1`)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]
pheatmap(decon.1.df[, ], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "T1D_MS_GSE60424 CIBERSORT meanGini 1220" )


idx <- grep("Whole", pdata$`celltype:ch1`)
decon.df_wholeblood <- decon.df[,idx]
pdata_wholeblood <- pdata[idx,]

mat_col <- data.frame(
                      Disease = pdata_wholeblood$`diseasestatus:ch1`
)
rownames(mat_col) <- colnames(decon.df_wholeblood)
mat_colors <- list(Disease = colorRampPalette( brewer.pal(8, "Set1"))(6)
)
names(mat_colors$Disease) <- unique( pdata_wholeblood$`diseasestatus:ch1`)
decon.1.df <- decon.df_wholeblood[rowMeans(decon.df_wholeblood)>=0.01,]
pheatmap(decon.1.df[, order(pdata_wholeblood$`diseasestatus:ch1`)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "T1D_MS_GSE60424 CIBERSORT meanGini 1220" )




idx <- grep("Whole", pdata$`celltype:ch1`, invert = T)
decon.df_wholeblood <- decon.df[,idx]
pdata_wholeblood <- pdata[idx,]

mat_col <- data.frame(CellType = pdata_wholeblood$`celltype:ch1`, 
                      Disease = pdata_wholeblood$`diseasestatus:ch1`
)
rownames(mat_col) <- colnames(decon.df_wholeblood)
mat_colors <- list(CellType = colorRampPalette( brewer.pal(6, "Accent"))(6), 
                   Disease = colorRampPalette( brewer.pal(6, "Set1"))(6)
)
names(mat_colors$CellType) <- unique( pdata_wholeblood$`celltype:ch1`)
names(mat_colors$Disease) <- unique( pdata_wholeblood$`diseasestatus:ch1`)
decon.1.df <- decon.df_wholeblood[rowMeans(decon.df_wholeblood)>=0.01,]
pheatmap(decon.1.df[, order(pdata_wholeblood$`celltype:ch1`)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "T1D_MS_GSE60424 CIBERSORT meanGini 1220" )

decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("matureNeutrophil|PMN", rownames(decon.1.df)),]) )
decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("Mono_C|Mono_I", rownames(decon.1.df)),]) )
decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("Bex|Bnaive|Bnsm", rownames(decon.1.df)),]) )
decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("CD4|Th1|Th2|Th17", rownames(decon.1.df)),]) )
decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("CD8|VD2", rownames(decon.1.df)),]) )
rownames(decon.1.df)[18:22] <- c("Neutrophils", "Monocytes", "B", "CD4T", "CD8T")
decon.1.df <- decon.1.df[18:22,]


#remove NK
idx <- grep("NK", pdata_wholeblood$`celltype:ch1`, invert = T)
decon.1.df <- decon.1.df[,idx]
pdata_wholeblood <- pdata_wholeblood[idx,]

#remove Sepsis 
idx <- grep("Sepsis", pdata_wholeblood$`diseasestatus:ch1`, invert = T)
decon.1.df <- decon.1.df[,idx]
pdata_wholeblood <- pdata_wholeblood[idx,]


mat_col <- data.frame(CellType = pdata_wholeblood$`celltype:ch1`, 
                      Disease = pdata_wholeblood$`diseasestatus:ch1`
)
rownames(mat_col) <- colnames(decon.1.df)
mat_colors <- list(CellType = colorRampPalette( brewer.pal(5, "Accent"))(5), 
                   Disease = colorRampPalette( brewer.pal(6, "Set1"))(5)
)
names(mat_colors$CellType) <- unique( pdata_wholeblood$`celltype:ch1`)
names(mat_colors$Disease) <- unique( pdata_wholeblood$`diseasestatus:ch1`)
pdata_wholeblood$`celltype:ch1` <- factor(pdata_wholeblood$`celltype:ch1`, levels = c("Neutrophils", "Monocytes", "B-cells", "CD4", "CD8"))
pheatmap(decon.1.df[, order(pdata_wholeblood$`celltype:ch1`)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = NA, 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "External Validation" )

melody.df <- decon.1.df
melody.df <- melt(melody.df)


##accuracy and RMSE
library(caret) 
library(Metrics)

idx <- grep("Neutrophils", pdata_wholeblood$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[1,idx])
colnames(acc.df)[2] <- "pred" 
#confusionMatrix(table(acc.df))
rmse.df <- rmse(acc.df$actual, acc.df$pred)

idx <- grep("Monocytes", pdata_wholeblood$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[2,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

idx <- grep("B-cells", pdata_wholeblood$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[3,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

idx <- grep("CD4", pdata_wholeblood$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[4,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

idx <- grep("CD8", pdata_wholeblood$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[5,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

rmse.df <- as.data.frame(rmse.df)
rmse.df <- rbind(rmse.df, colMeans(rmse.df))
rmse.df$cell <-  c("Neutrophils", "Monocytes", "B", "CD4T", "CD8T", "Average")
rmse.df$cell <- factor(rmse.df$cell, levels =  c("Average","Neutrophils", "Monocytes", "B", "CD4T", "CD8T" ) )
ggplot(rmse.df, aes(fill=cell, y=rmse.df, x=cell)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values = c("black", colorRampPalette( brewer.pal(5, "Accent"))(5))) +
  ggtitle("Validation Cohort") +ylab("RMSE")
rmse.df.benchmark <- rmse.df


###behcnmarking with CIBERSORT####
load(file = "./dataset/COVID/T1D_MS_GSE60424/T1D_MS_GSE60424.RData")
#rownames(tpm) <- gsub("-", "_", rownames(tpm))
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
sig_matrix <- read.delim(sig_matrix)
keep <- intersect(sig_matrix$Gene.symbol, rownames(tpm))
tpm.decon <- tpm[keep,]
pdata
idx <- grep("Whole", pdata$`celltype:ch1`, invert = T)
tpm.decon <- tpm.decon[,idx]
pdata <- pdata[idx,]
idx <- grep("NK", pdata$`celltype:ch1`, invert = T)
tpm.decon <- tpm.decon[,idx]
pdata <- pdata[idx,]
idx <- grep("Sepsis", pdata$`diseasestatus:ch1`, invert = T)
tpm.decon <- tpm.decon[,idx]
pdata <- pdata[idx,]

sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
decon.cibersort <- cibersort(sig_matrix, as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

#decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("Neutrophils", rownames(decon.1.df)),]) )
#decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("Mono_C|Mono_I", rownames(decon.1.df)),]) )
decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("B cells", rownames(decon.1.df)),]) )
decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("CD4|helper|Tregs", rownames(decon.1.df)),]) )
#decon.1.df <- rbind(decon.1.df, colSums(decon.1.df[grep("CD8", rownames(decon.1.df)),]) )
decon.1.df <- decon.1.df[c(10,9,11,12,3),]
rownames(decon.1.df) <- c("Neutrophils", "Monocytes", "B", "CD4T", "CD8T")

mat_col <- data.frame(CellType = pdata$`celltype:ch1`, 
                      Disease = pdata$`diseasestatus:ch1`
)
rownames(mat_col) <- colnames(decon.1.df)
mat_colors <- list(CellType = colorRampPalette( brewer.pal(5, "Accent"))(5), 
                   Disease = colorRampPalette( brewer.pal(6, "Set1"))(5)
)
names(mat_colors$CellType) <- unique( pdata$`celltype:ch1`)
names(mat_colors$Disease) <- unique( pdata$`diseasestatus:ch1`)
pdata$`celltype:ch1` <- factor(pdata$`celltype:ch1`, levels = c("Neutrophils", "Monocytes", "B-cells", "CD4", "CD8"))
pheatmap(decon.1.df[, order(pdata$`celltype:ch1`)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = F,show_rownames     = T, border_color = NA, 
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "External Validation" )

cibersort.df <- decon.1.df
cibersort.df <- melt(cibersort.df)
library(ggplot2)
library(ggpubr)
cor.df <- cbind(melody.df, cibersort.df)
cor.df <- cor.df[,c(1,3,6)]
colnames(cor.df) <- c("cell", "MELODY", "CIBERSORT")

ggscatter(cor.df, x = "MELODY", y = "CIBERSORT",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ 
  facet_wrap(~cell)+
  stat_cor(method = "pearson", label.x = 40, label.y = 40, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "Flow (Percentage)") + xlim(5,70) + ylim(25,100)



###RMSE
idx <- grep("Neutrophils", pdata$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[1,idx])
colnames(acc.df)[2] <- "pred" 
#confusionMatrix(table(acc.df))
rmse.df <- rmse(acc.df$actual, acc.df$pred)

idx <- grep("Monocytes", pdata$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[2,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

idx <- grep("B-cells", pdata$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[3,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

idx <- grep("CD4", pdata$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[4,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

idx <- grep("CD8", pdata$`celltype:ch1`)
acc.df <- data.frame(actual = rep(1, ncol(decon.1.df[,idx])))
acc.df <- cbind(acc.df, decon.1.df[5,idx])
colnames(acc.df)[2] <- "pred" 
rmse.df <- c(rmse.df,rmse(acc.df$actual, acc.df$pred) )

rmse.df <- as.data.frame(rmse.df)
rmse.df <- rbind(rmse.df, colMeans(rmse.df))
rmse.df$cell <-  c("Neutrophils", "Monocytes", "B", "CD4T", "CD8T", "Average")
rmse.df$cell <- factor(rmse.df$cell, levels =  c("Average","Neutrophils", "Monocytes", "B", "CD4T", "CD8T" ) )
ggplot(rmse.df, aes(fill=cell, y=rmse.df, x=cell)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values = c("black", colorRampPalette( brewer.pal(5, "Accent"))(5))) +
  ggtitle("Validation Cohort") +ylab("RMSE")

rmse.df.benchmark <- cbind(rmse.df.benchmark, rmse.df)
rmse.df.benchmark <- rmse.df.benchmark[,1:3]
colnames(rmse.df.benchmark) <- c("MELODY", "cell", "CIBERSORT")

ggplot(rmse.df.benchmark, aes(fill=cell, y=rmse.df, x=cell)) + 
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual(values = c("black", colorRampPalette( brewer.pal(5, "Accent"))(5))) +
  ggtitle("Validation Cohort") +ylab("RMSE")


##AVSc_GSE218474 ####
load(file = "./dataset/COVID/AVSc_GSE218474/AVSc_GSE218474.RData")

library(CIBERSORT)

rownames(tpm) <- gsub("-", "_", rownames(tpm))
pdata

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

colnames(pdata)[3:6] <- c("acs", "age", "avsc", "sex")
mat_col <- data.frame(ACS = pdata$acs, 
                      AVSc = pdata$avsc, 
                      Gender = pdata$sex
)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(ACS = colorRampPalette( brewer.pal(8, "Accent"))(2), 
                   AVSc = colorRampPalette( brewer.pal(8, "Set2"))(2), 
                   Gender = colorRampPalette( brewer.pal(8, "Set1"))(2)
)
names(mat_colors$ACS) <- unique( pdata$acs)
names(mat_colors$AVSc) <- unique( pdata$avsc)
names(mat_colors$Gender) <- unique( pdata$sex)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.1.df[, order(pdata$avsc)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "AVSc_GSE218474 CIBERSORT meanGini 1220" )

decon.1.df_avsc <- decon.1.df[,grep("Yes", pdata$avsc)]
pheatmap(decon.1.df_avsc[, order(pdata$acs[grep("Yes", pdata$avsc)])], 
         cluster_cols = T,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "T1D_MS_GSE60424 CIBERSORT meanGini 1220" )

decon.1.df_avsc.no <- decon.1.df[,grep("No", pdata$avsc)]
pheatmap(decon.1.df_avsc.no[, order(pdata$acs[grep("No", pdata$avsc)])], 
         cluster_cols = T,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "T1D_MS_GSE60424 CIBERSORT meanGini 1220" )

pdata <- cbind(pdata, t(decon.1.df))
library(ggplot2)
library(ggpubr)
pdata%>%
  ggplot( aes(x=avsc, y=PMN, fill=avsc)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("No", "Yes")), y_position = 0.5)

pdata%>%
  ggplot( aes(x=acs, y=PMN, fill=acs)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("avsc")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("NSTEMI", "STEMI")), y_position = 0.5)

pdata%>%
  ggplot( aes(x=avsc, y=Mono_C, fill=avsc)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_C ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("No", "Yes")), y_position = 0.3)
pdata%>%
  ggplot( aes(x=avsc, y=Mono_I, fill=avsc)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_I ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("No", "Yes")), y_position = 0.3)

pdata%>%
  ggplot( aes(x=avsc, y=Tfh, fill=avsc)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Tfh ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("No", "Yes")), y_position = 0.3)







##AMI_GSE249812 ####
load(file = "./dataset/COVID/AMI_GSE249812/AMI_GSE249812.RData")
sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]
library(CIBERSORT)

rownames(tpm) <- gsub("-", "_", rownames(tpm))
pdata

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

colnames(pdata)[3] <- c("group")
mat_col <- data.frame(group = pdata$group
)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(group = colorRampPalette( brewer.pal(8, "Accent"))(3)
                   )
names(mat_colors$group) <- unique( pdata$group)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.1.df[, ], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "AMI_GSE249812 CIBERSORT meanGini 1220" )


pdata <- cbind(pdata, t(decon.1.df))
pdata$group <- factor(pdata$group, levels = c("acute stage", "subacute stage", "old phase of AMI"))
library(ggplot2)
library(ggpubr)
pdata%>%
  ggplot( aes(x=group, y=CD8Teff, fill=group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=5, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F", "darkred"))+
  #facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("acute stage", "subacute stage")), y_position = 0.5)
pdata%>%
  ggplot( aes(x=group, y=Th2, fill=group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=5, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F", "darkred"))+
  #facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("acute stage", "subacute stage")), y_position = 0.3)






##AMI_GSE166780 ####
load(file = "./dataset/COVID/AMI_GSE166780/AMI_GSE166780.RData")

sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]
library(CIBERSORT)

rownames(tpm) <- gsub("-", "_", rownames(tpm))

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

colnames(pdata)[4] <- c("group")
mat_col <- data.frame(group = pdata$group
)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(group = colorRampPalette( brewer.pal(8, "Accent"))(3)
)
names(mat_colors$group) <- unique( pdata$group)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.1.df[, ], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "AMI_GSE249812 CIBERSORT meanGini 1220" )


##STEMI_GSE103182####
load(file = "./dataset/COVID/STEMI_GSE103182/STEMI_GSE103182.RData")
#sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]
library(CIBERSORT)

rownames(tpm) <- gsub("-", "_", rownames(tpm))

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

colnames(pdata)[3] <- c("group")
mat_col <- data.frame(group = pdata$group
)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(group = colorRampPalette( brewer.pal(3, "Accent"))(2)
)
names(mat_colors$group) <- unique( pdata$group)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.1.df[, ], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "STEMI_GSE103182 CIBERSORT meanGini 1220" )


pdata <- cbind(pdata, t(decon.1.df))
#pdata$group <- factor(pdata$group, levels = c())
library(ggplot2)
library(ggpubr)
p1 <- pdata%>%
  ggplot( aes(x=group, y=CD8Tcm, fill=group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=3, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("CD8Tcm ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  #facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("STEMI", "NSTEMI")), y_position = 0.15)
p2 <- pdata%>%
  ggplot( aes(x=group, y=MAIT, fill=group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=3, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("MAIT ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  #facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("STEMI", "NSTEMI")), y_position = 0.06)
p3 <- pdata%>%
  ggplot( aes(x=group, y=VD2, fill=group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=3, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("VD2 ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  #facet_wrap("acs")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("STEMI", "NSTEMI")), y_position = 0.15)
ggarrange(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)  




##Spesis_PRJCA006118####
load(file = "./dataset/COVID/Spesis_PRJCA006118/Spesis_PRJCA006118.RData")
survival <- read.csv("./dataset/COVID/Spesis_PRJCA006118/OMIX004883-01.csv", row.names = 1)
#sig.mat <- sig.mat[,grep("matureNeutrophil|Erythroblast|Megakaryocyte" , colnames(sig.mat), invert = T)]
pdata$id <- gsub("_SRM", "", pdata$id)
tmp <- data.frame(matrix(unlist(strsplit(pdata$id, split = "_")), ncol =  3, byrow = T ) )
colnames(tmp) <- c("center", "pID", "day")
pdata <- cbind(pdata, tmp)
pdata$Treatment <- gsub(0, "No", pdata$Treatment)
pdata$Treatment <- gsub("1", "Yes", pdata$Treatment)
rownames(pdata) <- pdata$id
pdata$mortality <- NA

mort_flg <- survival[which(rownames(survival) %in% rownames(pdata)),"mort_flg"]
pdata$mortality[which(rownames(pdata) %in% rownames(survival))] <- mort_flg

library(CIBERSORT)

rownames(tpm) <- gsub("-", "_", rownames(tpm))

library(limma)
plotMDS(log2(tpm+1), top = 2000, gene.selection = "common", pch = 19 )

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]), QN = F)
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])


mat_col <- data.frame(treatment = pdata$Treatment, 
                      #center = pdata$center, 
                      day = pdata$day, 
                      mortality = pdata$mortality,
                      SOFA = pdata$SOFA
)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(treatment = colorRampPalette( brewer.pal(3, "Accent"))(2), 
                   day = colorRampPalette( brewer.pal(3, "Dark2"))(3),
                   mortality = colorRampPalette( brewer.pal(3, "Set1"))(3)[c(3,2,1)],
                   SOFA = colorRampPalette( brewer.pal(9, "YlOrRd"))(21)
)
names(mat_colors$treatment) <- unique( pdata$Treatment)
names(mat_colors$day) <- unique( pdata$day)
names(mat_colors$SOFA) <- unique( pdata$SOFA)
names(mat_colors$mortality) <- unique( pdata$mortality)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]

pheatmap(decon.df[, order(pdata$SOFA) ], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Spesis_PRJCA006118 CIBERSORT meanGini 1220" )

pheatmap(decon.df[, order(pdata$mortality) ], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Spesis_PRJCA006118 CIBERSORT meanGini 1220" )

idx <- grep("d1", pdata$day)
pdata.d1 <- pdata[idx,]
decon.d1 <- decon.df[,idx]
pheatmap(decon.d1[, order(pdata.d1$mortality) ], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Spesis_PRJCA006118 CIBERSORT meanGini 1220" )


idx <- grep("d3", pdata$day)
pdata.d1 <- pdata[idx,]
decon.d1 <- decon.df[,idx]
pheatmap(decon.d1[, order(pdata.d1$mortality) ], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Spesis_PRJCA006118 CIBERSORT meanGini 1220" )

