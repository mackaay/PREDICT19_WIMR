library(dplyr)
library(Seurat)
#library(patchwork)
library(RColorBrewer)

pbmc.data <- Read10X(data.dir = "./dataset/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc.sce <- as.SingleCellExperiment(pbmc)

#Cell type annotation 
library(celldex)
ref <- BlueprintEncodeData()
ref
#ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

library(SingleR)
pred <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.main, assay.type.test=1)
table(pred$labels)
pbmc.sce$cell.main <- pred$labels
pred <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred$labels)
pbmc.sce$cell.fine <- pred$labels


###testing seurat PBMC data####
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/RF_meanGini.rds")
niteration <- 1220
sig.mat_pbmc <- sig.mat[,grep("Basophils|Erythroblast|matureNeutrophil", colnames(sig.mat), invert = T)]

ref.bulk <- log2(sig.mat[,]+1)
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred.bulk$labels)
pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)

pred.df <- data.frame(SingleR = pred.singler$labels, bulk = pred.bulk$labels, 
                      freq = c(1))
pred.df$bulk <- factor(pred.df$bulk, levels = c("CD4Tnaive", "CD4Teff", "Tfh", "Th1", "Th2", "Th17", 
                                                "CD8Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem", "MAIT", "VD2", 
                                                "Mono_C" , "Mono_I", "Mono_NC", "Bnaive", "Bex", "Bnsm", "Bsm", 
                                                "NK", "Plasmablasts", "mDC", "pDC","Treg", "PMN"))
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
  theme_void() + ggtitle("seurat pbmc dataset") +  theme(legend.position="none")
sankeyp



####clustering#####
library(ggplot2)
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
DimPlot(pbmc, label = TRUE, reduction = "umap") + NoLegend() +ggtitle("Seurat clustering")

pbmc$singler.fine <- pred.singler$labels
pbmc$bulk <- pred.bulk$labels


pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(17) 
Idents(pbmc) <- 'RNA_snn_res.0.3'
p2 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()+ggtitle("Seurat clustering")
pal.celltype <- brewer.pal(8, "Set2")
pal.celltype <- colorRampPalette(pal.celltype)(38) 
Idents(pbmc) <- 'singler.fine'
p3 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend() +ggtitle("SingleR")
pal.celltype <- brewer.pal(8, "Set3")
pal.celltype <- colorRampPalette(pal.celltype)(26) 
Idents(pbmc) <- 'bulk'
p4 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()+ggtitle("MELODY")
p2 + p3 + p4


FeaturePlot(pbmc, features = c("PPBP", "LY6G", "LY6C", "CD11B", "CD45", "CD14"))


###testing Jin's PMN data####
load("/datasets/work/hb-diab-cfdna/work/scratch/chenkai/scRNA/scRNA_MDSC_mm10/pmn_rename.rds")
pmn.sce <- as.SingleCellExperiment(pmn)
rownames(pmn.sce) <- toupper(rownames(pmn.sce))

pred <- SingleR(test=pmn.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred$labels)

pred <- SingleR(test=pmn.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred$labels)

Idents(pmn)




###COVID and IAV ####
#https://www.cell.com/immunity/pdf/S1074-7613(20)30316-2.pdf
#COVID and IAV
pbmc <- readRDS("./dataset/COVID/Final_nCoV_0716_upload.RDS")
pbmc.sce <- as.SingleCellExperiment(pbmc)

pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred.singler$labels)
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred.bulk$labels)

pbmc$singler.fine <- pred.singler$labels
pbmc$bulk <- pred.bulk$labels

pal.celltype <- brewer.pal(8, "Set1")

pal.celltype <- colorRampPalette(pal.celltype)(15) 
Idents(pbmc) <- 'cell_type'
p1 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- colorRampPalette(pal.celltype)(19) 
Idents(pbmc) <- 'integrated_snn_res.0.8'
p2 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- colorRampPalette(pal.celltype)(22) 
Idents(pbmc) <- 'singler.fine'
p3 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- colorRampPalette(pal.celltype)(29) 
Idents(pbmc) <- 'bulk'
p4 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
p1 + p2 + p3 + p4

Idents(pbmc) <- 'integrated_snn_res.0.8'
tcell <- pbmc[,grep("0|1|2|3|4|5|9|10|16", Idents(pbmc))]

pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(15) 
Idents(tcell) <- 'cell_type'
p1 <- DimPlot(tcell, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(19) 
Idents(tcell) <- 'integrated_snn_res.0.8'
p2 <- DimPlot(tcell, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(22) 
Idents(tcell) <- 'singler.fine'
p3 <- DimPlot(tcell, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(26) 
Idents(tcell) <- 'bulk'
p4 <- DimPlot(tcell, cols = pal.celltype, label = T) + NoLegend()
p1 + p2 + p3 + p4


###severeCOVIDFlu_GSE149689 ####
pbmc.data <- Read10X(data.dir = "./dataset/COVID/severeCOVIDFlu_GSE149689/")
pbmc <- CreateSeuratObject(counts = pbmc.data,  min.cells = 3, min.features = 200)
pbmc
pbmc$group <- "nCoV"
pbmc$group[grep("-3$|-4$|-6$|-7$|-8$",colnames(pbmc))] <- "Flu"
pbmc$group[grep("-5$|-13$|-14$|-19$",colnames(pbmc))] <- "Flu"
pbmc.sce <- as.SingleCellExperiment(pbmc)

pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred.singler$labels)
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred.bulk$labels)

pbmc$singler.fine <- pred.singler$labels
pbmc$bulk <- pred.bulk$labels

pred.df <- data.frame(SingleR = pred.singler$labels, bulk = pred.bulk$labels, 
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



###COVID_284_GSE158055####
library(GEOquery)
GSEset <- getGEO("GSE158055", destdir = "./dataset/COVID/COVID_284_GSE158055/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 55:71)]
pdata <-pdata[grep("PBMC", pdata$`sample type:ch1`),]
rownames(pdata) <- pdata$title

pbmc <- Read10X(data.dir = "./dataset/COVID/COVID_284_GSE158055/part1/", gene.column = 1)
pbmc <- CreateSeuratObject(counts = pbmc, project = "COVID_284_part1", min.cells = 10, min.features = 200)
pbmc


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
DimPlot(pbmc, label = TRUE, reduction = "umap") + NoLegend() +ggtitle("Seurat clustering")

pbmc.sce <- as.SingleCellExperiment(pbmc)
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred.bulk$labels)
pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)

length(pred.singler$labels)
pbmc$singler.fine <- pred.singler$labels
pbmc$bulk <- pred.bulk$labels


pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(24) 
Idents(pbmc) <- 'RNA_snn_res.0.3'
p2 <- DimPlot(pbmc, cols = pal.celltype, label = T, raster = F) + NoLegend()+ggtitle("Seurat clustering")
pal.celltype <- brewer.pal(8, "Set2")
pal.celltype <- colorRampPalette(pal.celltype)(38) 
Idents(pbmc) <- 'singler.fine'
p3 <- DimPlot(pbmc, cols = pal.celltype, label = T, raster = F) + NoLegend() +ggtitle("SingleR")
pal.celltype <- brewer.pal(8, "Set3")
pal.celltype <- colorRampPalette(pal.celltype)(29) 
Idents(pbmc) <- 'bulk'
p4 <- DimPlot(pbmc, cols = pal.celltype, label = T, raster = F) + NoLegend()+ggtitle("MELODY")
p2 + p3 + p4

library(caret)
library(pheatmap)
conf.mat <- as.matrix(table(pbmc$singler.fine, pbmc$bulk))

pheatmap(as.matrix(table(pbmc$singler.fine, pbmc$bulk)), scale = "column",
         cluster_cols = F,cluster_rows = F, display_numbers = T, color = colorRampPalette(c( "white", "#FFFF37","#FF9700","#FF3A00", "#FF0000"))(100), 
         show_colnames     = T,show_rownames     = T, border_color = "black", 
         #annotation_col    = mat_col,annotation_colors = mat_colors,
         main = c("Confusion Matrix" ))







###GAC cancer and LN and PBMC and etc ####
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239676
pbmc.data <- Read10X(data.dir = "./dataset/COVID/GAC_GSE239676/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GAC", min.cells = 3, min.features = 200)
pbmc



###ESRD PBMC ####
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233315
pbmc.data <- Read10X(data.dir = "./dataset/COVID/ESRD_GSE233315/ESRD//")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ESRD", min.cells = 3, min.features = 200)
pbmc
pbmc.data <- Read10X(data.dir = "./dataset/COVID/ESRD_GSE233315/NC//")
tmp <- CreateSeuratObject(counts = pbmc.data, project = "NC", min.cells = 3, min.features = 200)
keep <- intersect(rownames(pbmc), rownames(tmp))
pbmc <- merge(pbmc[keep,], y = tmp[keep,])

pbmc.sce <- as.SingleCellExperiment(pbmc)

pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred.singler$labels)
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred.bulk$labels)
system.time(SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))) ; system.time(SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1))

pbmc$singler.fine <- pred.singler$labels
pbmc$bulk <- pred.bulk$labels

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
pal.celltype <- colorRampPalette(pal.celltype)(17) 
Idents(pbmc) <- 'RNA_snn_res.0.3'
p2 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(38) 
Idents(pbmc) <- 'singler.fine'
p3 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(26) 
Idents(pbmc) <- 'bulk'
p4 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
p2 + p3 + p4


###Bladder Cancer PBMC####
#https://rupress.org/jem/article/221/8/e20240045/276792/Urine-scRNAseq-reveals-new-insights-into-the
pbmc.data <- Read10X(data.dir = "/datasets/work/lw-project-airway/work/User/chenkai/exvivo_RNAseq/dataset/COVID/BladderCancer_GSE267718/Patient5PBMC")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "Patient6", min.cells = 3, min.features = 200)
pbmc
pbmc.sce <- as.SingleCellExperiment(pbmc)

pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred.singler$labels)
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred.bulk$labels)

pbmc$singler.fine <- pred.singler$labels
pbmc$bulk <- pred.bulk$labels

pred.df <- data.frame(SingleR = pred.singler$labels, bulk = pred.bulk$labels, 
                      freq = c(1))
#sankey plot
library(ggalluvial)
ggplot(data = pred.df,
       aes(axis1 = SingleR, axis2 = bulk, y = freq)) +
  geom_alluvium(aes(fill = bulk)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("SingleR", "bulk"),
                   expand = c(0.15, 0.05)) +
  theme_void()

pbmc$singler.fine <- pred.singler$labels
pbmc$bulk <- pred.bulk$labels
table(pbmc$cell_type)


pbmc <- pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() #features = rownames(pbmc.combined)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims = 30)
pbmc <- FindNeighbors(object = pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = c( 0.1, 0.3))
pbmc <- RunUMAP(pbmc, reduction = "pca",  dims = 1:20)

pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(11) 
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
p2 + p3 + p4






###MIA scRNAseq  FFPE tissue####
#https://www.cell.com/cell-reports/fulltext/S2211-1247(24)00720-4
library(ggplot2)
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

pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred.singler$labels)
pbmc$singler.fine <- pred.singler$labels
pbmc <- pbmc[,grep("Melanocytes|Endo", pbmc$singler.fine, invert = T)] #remove melanoma cells
pbmc.sce <- as.SingleCellExperiment(pbmc)

pred.singler <- SingleR(test=pbmc.sce, ref=ref, labels=ref$label.fine, assay.type.test=1)
table(pred.singler$labels)
pbmc$singler.fine <- pred.singler$labels
pred.bulk <- SingleR(test=pbmc.sce, ref=ref.bulk, labels=colnames(ref.bulk))
table(pred.bulk$labels)
pbmc$bulk <- pred.bulk$labels
pred.df <- data.frame(SingleR = pred.singler$labels, bulk = pred.bulk$labels, 
                      freq = c(1))

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
DotPlot(pbmc, features = c("EPCAM", "CDH1"))

pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(15) 
Idents(pbmc) <- 'RNA_snn_res.0.3'
p2 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(38) 
Idents(pbmc) <- 'singler.fine'
p3 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(26) 
Idents(pbmc) <- 'bulk'
p4 <- DimPlot(pbmc, cols = pal.celltype, label = T) + NoLegend()
p2 + p3 + p4







