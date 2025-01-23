library(dplyr)
library(Seurat)
#library(patchwork)
library(RColorBrewer)
library(Matrix)
#https://www.sciencedirect.com/science/article/pii/S0092867421001483?via%3Dihub#sec6



library(GEOquery)
GSEset <- getGEO("GSE158055", destdir = "./dataset/COVID/COVID_284_GSE158055/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 55:71)]
pdata <-pdata[grep("PBMC", pdata$`sample type:ch1`),]
rownames(pdata) <- pdata$title

cell.anno <- read.csv("./dataset/COVID/COVID_284_GSE158055/GSE158055_cell_annotation.csv.gz")
table(cell.anno$majorType, cell.anno$celltype)
cell.anno <- cell.anno[grep("B|CD|DC|Mega|Mono|Neu|NK|Plasma", cell.anno$majorType),]
unique(cell.anno$sampleID)[1]

#part 1
pbmc.data <- Read10X(data.dir = "./dataset/COVID/COVID_284_GSE158055/part1/", gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "COVID_284_part1", min.cells = 10, min.features = 200)
pbmc

keep <- intersect(colnames(pbmc), cell.anno$cellName)
pbmc <- pbmc[,keep]
rownames(cell.anno) <- cell.anno$cellName
cell.anno_1 <- cell.anno[keep,]
pbmc$sampleID <- cell.anno_1$sampleID
pbmc$majorType <- cell.anno_1$majorType
actual <- as.data.frame(table(pbmc$sampleID, pbmc$majorType))

pbmc <- pbmc %>%
  NormalizeData()
pseudobulk <- AggregateExpression(pbmc, assays = "RNA", return.seurat = F, group.by = "sampleID")
pseudobulk <- pseudobulk$RNA

keep <- intersect(colnames(pseudobulk), rownames(pdata))
pseudobulk_1 <- pseudobulk[,keep]
pdata_1 <- pdata[keep,]



#part 2
pbmc.data <- Read10X(data.dir = "./dataset/COVID/COVID_284_GSE158055/part2/", gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "COVID_284_part2", min.cells = 10, min.features = 200)
pbmc
#keep <- intersect(colnames(pbmc), colnames(pbmc2))
#pbmc2 <- pbmc2[,!colnames(pbmc2) %in% keep]
#pbmc <- merge(pbmc, y = pbmc2, project = "COVID_284")
keep <- intersect(colnames(pbmc), cell.anno$cellName)
pbmc <- pbmc[,keep]
rownames(cell.anno) <- cell.anno$cellName
cell.anno_2 <- cell.anno[keep,]
pbmc$sampleID <- cell.anno_2$sampleID
pbmc$majorType <- cell.anno_2$majorType
actual <- as.data.frame(table(pbmc$sampleID, pbmc$majorType))

pbmc <- pbmc %>%
  NormalizeData()
pseudobulk <- AggregateExpression(pbmc, assays = "RNA", return.seurat = F, group.by = "sampleID")
pseudobulk <- pseudobulk$RNA

keep <- intersect(colnames(pseudobulk), rownames(pdata))
pseudobulk_2 <- pseudobulk[,keep]
pdata_2 <- pdata[keep,]


cell.anno.merge <- rbind(cell.anno_1, cell.anno_2)
keep <- intersect(rownames(pseudobulk_1), rownames(pseudobulk_2))
pseudobulk.merge <- cbind(pseudobulk_1[keep,], pseudobulk_2[keep,])
pdata.merge <- rbind(pdata_1, pdata_2)


save(cell.anno.merge, pseudobulk.merge, pdata.merge, file = "./dataset/COVID/COVID_284_GSE158055/COVID_284_GSE158055.RData")


##CIBERSORT  ####
load(file = "./dataset/COVID/COVID_284_GSE158055/COVID_284_GSE158055.RData")
library(pheatmap)
library(RColorBrewer)
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/RF_meanGini.rds")
load(file = "./dataset/COVID/RF_meanGini2.rds")
load(file = "./dataset/COVID/gene_length.rds")
niteration <- 1220
library(CIBERSORT)

###convert to TPM####
library(edgeR)
y <- DGEList(counts = pseudobulk.merge)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 4)

keep <- intersect(rownames(all.tmm), names(gene.length))
all.tmm <- all.tmm[keep,]
gene.length <- gene.length[keep]
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)
rownames(tpm) <- gsub("-", "_", rownames(tpm))
rownames(pseudobulk.merge) <- gsub("-", "_", rownames(pseudobulk.merge))

colnames(sig.mat)
sig.mat_pbmc <- sig.mat[,grep("Basophils|Erythroblast|matureNeutrophil|Megakaryocyte", colnames(sig.mat), invert = T)]

keep <- intersect(rownames(sig.mat_pbmc[rownames(meanGini)[1:niteration],]), rownames(tpm))
tpm.decon <- pseudobulk.merge[keep,]
decon.cibersort <- cibersort(sig.mat_pbmc[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]))
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.cibersort
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

keep <- intersect(colnames(decon.df), rownames(pdata.merge))
pdata.merge <- pdata.merge[keep, ]
decon.df <- decon.df[,keep]


pdata.merge$`outcome:ch1`
mat_col <- data.frame(Severity = pdata.merge$`covid-19 severity:ch1`, 
                      Outcome = pdata.merge$`outcome:ch1`, 
                      Time = pdata.merge$`sample time:ch1`, 
                      Type = pdata.merge$`sample type:ch1`)
mat_col$Type <- gsub(" sorted from fresh PBMC (FACS)", "", mat_col$Type )
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Outcome = brewer.pal(3, "Set1")[1:3], 
                   Severity = brewer.pal(3, "YlOrRd")  , 
                   Time = brewer.pal(3, "Dark2"), 
                   Type = brewer.pal(7, "Accent"))
names(mat_colors$Outcome) <- unique(pdata.merge$`outcome:ch1`)
names(mat_colors$Severity) <-unique(pdata.merge$`covid-19 severity:ch1`)
names(mat_colors$Time) <- unique(pdata.merge$`sample time:ch1`)
names(mat_colors$Type) <-unique(pdata.merge$`sample type:ch1`)
#decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]
pheatmap(decon.df[,order(pdata.merge$`covid-19 severity:ch1`, decreasing = F)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID 248 cohort CIBERSORT meanGini 1220" )

table(cell.anno.merge$majorType)

decon.bcell <- decon.df[,grep("B cell", pdata.merge$`sample type:ch1`)]
pheatmap(decon.bcell[,], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID 248 cohort CIBERSORT meanGini 1220" )







library(ggplot2)
library(ggpubr)
actual[grep("S-HC001", actual$Var1),]
actual$Var1 <- as.character(actual$Var1)
actual$Var2 <- factor(actual$Var2, levels = c("Mega", "B", "Plasma", "CD4", "CD8", "NK", "DC", "Mono" , "Neu"))
p1 <- ggplot(actual[grep("S-HC001", actual$Var1),], aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") 
cibersort.pred  <- data.frame(Var1 = "S-HC001", Var2 = colnames(decon.cibersort)[1:(ncol(decon.cibersort)-3)], Freq = decon.cibersort["S-HC001", 1:(ncol(decon.cibersort)-3)]) 
cibersort.pred$Var2  <- factor(cibersort.pred$Var2, levels = c("Megakaryocyte", "Bnaive", "Bnsm", "Bsm", "Bex", "Plasmablasts", "CD4Tnaive", "CD4Teff", "Tfh", "Th1", "Th2", "Th17", "Treg",
                                                               "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem", "MAIT", "VD2", "NK", "mDC", "pDC", "Mono_C", "Mono_I", "Mono_NC", "PMN")) 
p2 <- ggplot(cibersort.pred, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") 
ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
