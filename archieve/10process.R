library(RColorBrewer)
library(edgeR)

counts <- read.delim("./dataset/COVID/gene_id_exon.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))

#FGL2
counts[grep("ENSG00000127951", rownames(counts)),] 
#RN7SL1
counts[grep("ENSG00000276168", rownames(counts)),] 


###mapping gene name ####
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
gene.length <- counts$Length
names(gene.length) <- rownames(counts)
counts <- counts[,-1]

save(gene.length, file = "./dataset/COVID/gene_length.rds")

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

boxplot(log10(tpm[,1:10]) )

colnames(tpm) <- gsub("_Whole_blood", "", colnames(tpm))
colnames(tpm) <- gsub("COVID_positive", "COVID", colnames(tpm))
colnames(tpm) <- gsub("healthy_control", "healthy", colnames(tpm))
pData <- matrix(as.character( unlist(strsplit(colnames(tpm), split = "_")) ), ncol=5, byrow = T)
pData <- as.data.frame(pData)
pData <- pData[,c(1,2,3,5)]

colnames(tpm) <- pData$V1

write.csv(tpm, file = "./tpm_COVIDcohort.csv")
write.table(tpm, file = "./tpm_COVIDcohort.txt", sep = "\t")
write.csv(tpm[,1:100], file = "./tpm_COVIDcohort1.csv"); write.csv(tpm[,101:200], file = "./tpm_COVIDcohort2.csv"); write.csv(tpm[,201:300], file = "./tpm_COVIDcohort3.csv"); write.csv(tpm[,301:400], file = "./tpm_COVIDcohort4.csv"); write.csv(tpm[,401:ncol(tpm)], file = "./tpm_COVIDcohort5.csv");
save(tpm, counts, rpkm, pData, file = "./COVIDcohort.RData")




###PBMC GEO ABIS####
library(RColorBrewer)
library(edgeR)

counts <- read.delim("./dataset/COVID/gene_id_exon_BP_pbmc_GEO_ABIS.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(2,3,4,5)]
rownames(pdata) <- pdata$V2
write.csv(pdata, file = "./tmp.csv")
pdata <- read.csv("./dataset/COVID/PBMC_GEO_pdata.csv", row.names = 1)


#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$id
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

boxplot(log10(tpm[,1:50]) )
boxplot(log10(rpkm[,1:50]) )


save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/PBMC_GEO_ABIS.RData")



###GMDSC_GSE178824 with COVID and control##### 
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_GMDSC_GSE178824.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(2,3,4,5)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$V2
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)


save(tpm, counts, rpkm, gene.length, pdata, file = "./dataset/COVID/PBMC_GEO_GMDSC_GSE178824.RData")



###M-MDSC_GSE199286 with COVID and control##### 
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_MMDSC_GSE199286.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
tmp$V4[grep("\\.", tmp$V4 , invert = T)] <- "neg"
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(1,2,3,4)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$V2
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)


save(tpm, counts, rpkm, gene.length, pdata, file = "./dataset/COVID/PBMC_GEO_MMDSC_GSE199286.RData")




###monoMMDSC_GSE183854##### 
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_monoMMDSC_GSE183854.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(1,2,3,4)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$V2
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)


save(tpm, counts, rpkm, gene.length, pdata, file = "./dataset/COVID/PBMC_GEO_monoMMDSC_GSE183854.RData")




###neutrophil_GSE163533##### 
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_neutrophil_GSE163533.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(1,2,3,4)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$V2
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)


save(tpm, counts, rpkm, gene.length, pdata, file = "./dataset/COVID/PBMC_GEO_neutrophil_GSE163533.RData")





###neutrophil_GSE66895##### 
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_neutrophil_GSE66895.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(1,2,3,4)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$V2
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)


save(tpm, counts, rpkm, gene.length, pdata, file = "./dataset/COVID/PBMC_GEO_neutrophil_GSE66895.RData")




###validation ##### 
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_validation.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(1,2,3,4)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$V2
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)
x <- counts / gene.length
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
tpm.mat <- round(tpm.mat, 4)

data.frame(tpm$GSM5970292, tpm.mat[,1])

save(tpm, counts, rpkm, gene.length, pdata, file = "./dataset/COVID/PBMC_GEO_validation.RData")




###PMN_GSE163834##### 
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_PMN_GSE163834.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(1,2,3,7)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
gene.length <- counts$Length
colnames(counts)[2:ncol(counts)] <- pdata$V2
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)


save(tpm, counts, rpkm, gene.length, pdata, file = "./dataset/COVID/PBMC_GEO_PMN_GSE163834.RData")



###Erythrocytes and Megakaryocytes####
library(limma)
library(celldex)
ref <- BlueprintEncodeData()
ref
ref[,grep("Erythrocyte|Megakaryocyte|Neutrophil",colData(ref)$label.fine)]
assays(ref[,grep("Erythrocyte|Megakaryocyte|Neutrophil",colData(ref)$label.fine)])[["logcounts"]]
Erythro_Mega_Neutro <- 2^assays(ref[,grep("Erythrocyte|Megakaryocyte|Neutrophil",colData(ref)$label.fine)])[["logcounts"]] -1
Erythro_Mega_Neutro <- Erythro_Mega_Neutro[,order(colnames(Erythro_Mega_Neutro))]
colnames(Erythro_Mega_Neutro)
pdata <- data.frame(id = colnames(Erythro_Mega_Neutro), cell = c(rep("bandNeutrophil", 3),  rep("Megakaryocyte", 5 ), rep("Erythroblast", 7), rep("matureNeutrophil", 17), rep("Neutrophil_BM",3)))

plotMDS(log2(Erythro_Mega_Neutro+1), top=1000, gene.selection="common", labels = pdata$cell)

#sample distance cluster 
cell_dist <- dist(t(log2(Erythro_Mega_Neutro[grep("MIR[0-9]|^MT-|^RPL[0-9]", rownames(Erythro_Mega_Neutro), invert = T),]+1)))
cell_hclust <- hclust(cell_dist, method = "ward.D")
# The default `plot()` function can be used to produce a simple dendrogram
plot(cell_hclust, label = pdata$cell, type = "triangle", ylab = "Height")

idx <- grep("BM", pdata$cell, invert = T)
Erythro_Mega_Neutro <- Erythro_Mega_Neutro[,idx]
pdata <- pdata[idx,]

cell_dist <- dist(t(log2(Erythro_Mega_Neutro[grep("MIR[0-9]|^MT-|^RPL[0-9]", rownames(Erythro_Mega_Neutro), invert = T),]+1)))
cell_hclust <- hclust(cell_dist, method = "ward.D")
# The default `plot()` function can be used to produce a simple dendrogram
plot(cell_hclust, label = pdata$cell, type = "triangle", ylab = "Height")

save(Erythro_Mega_Neutro, pdata, file = "./dataset/COVID/Erythro_Mega_Neutro.RData")

####
mat <- matrix(1:10, nrow = 2)
sweep(mat, 1, rowSums(mat), `/`)
mat/ colSums(mat)
####




###GTEx dataset ####
gtex <- read.delim(file="/datasets/work/lw-project-airway/work/Data/level2/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", skip=2)



###Japan COVID cohort GSE206263####
counts <- read.delim("./dataset/COVID/gene_id_exon_BP_COVID_GSE206263.txt", stringsAsFactors = F, skip = 1)
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))
tmp <- counts[,grep("GSM6248833", colnames(counts))]
counts <- counts[, grep("GSM6248833", colnames(counts), invert = T)]
counts$SRR19681197_GSM6248833_180E1A_14_Severe_Positive <-  rowSums(tmp)
tmp <- as.data.frame(sapply( strsplit(colnames(counts)[2:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
colnames(counts[2:ncol(counts)]) <- tmp$V2
pdata <- tmp[,c(2,3,4,5,6)]
rownames(pdata) <- pdata$V2

#mapping gene name
library(EnsDb.Hsapiens.v86)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
gene.length <- counts$Length
colnames(pdata) <- c("id", "sample", "days", "group", "test")
colnames(counts)[2:ncol(counts)] <- pdata$id
counts <- counts[,-1]

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
tmp <- sweep(counts, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

boxplot(log10(tpm[,1:50]) )
boxplot(log10(rpkm[,1:50]) )


save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/PBMC_JapanCOVID_GSE206263.RData")


###COVID cohort GSE215865####
counts <- read.csv("./dataset/COVID/GSE215865_rnaseq_raw_count_matrix.csv.gz")
#counts <- read.delim("./dataset/gene_id1.txt", stringsAsFactors = F, skip = 1)
counts$Ensembl_Gene_ID <- substr(counts$Ensembl_Gene_ID, start = 1L, stop = 15L)
counts <- aggregate(x = counts[,2:ncol(counts)],by = list(counts$Ensembl_Gene_ID), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

tmp <- as.data.frame(sapply( strsplit(colnames(counts)[1:ncol(counts)], split = "_"), `[` , 1:10 ) )
tmp <- as.data.frame(t(tmp))
#colnames(counts[1:ncol(counts)]) <- tmp$V2

library(GEOquery)
GSEset <- getGEO("GSE215865", destdir = "./dataset/COVID/COVID_GSE215865/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 62:74)]
rownames(pdata) <- paste(pdata$`blood sample_id:ch1`, pdata$`library prep_plate:ch1`, sep = "_")

keep <- intersect(rownames(pdata), colnames(counts))
pdata <- pdata[keep, ]
counts <- counts[, keep]


#mapping gene name
library(EnsDb.Hsapiens.v86)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
             )
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load(file = "./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))

all.tmm <- all.tmm[keep, ]
data.frame(rownames(all.tmm), names(gene.length))
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/WBC_COVID_GSE215865.RData")


###HBV GSE173897#####
counts <- read.csv("./dataset/COVID/HBV_GSE173897/gene_id_exon_BP_HBV.txt", stringsAsFactors = F, skip = 1, sep = "\t")
rownames(counts) <- counts$Geneid
counts <- counts[,6:ncol(counts)]
colnames(counts) <- gsub("\\..alignment.", "", colnames(counts))
colnames(counts) <- gsub("_s.bam", "", colnames(counts))

#mapping gene name
library(EnsDb.Hsapiens.v86)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(counts), keytype = "GENEID", columns = c("SYMBOL","GENEID", "GENEBIOTYPE"))
rownames(geneIDs) <- geneIDs$GENEID
keep <- intersect(rownames(geneIDs), rownames(counts))
counts <- counts[keep,]
geneIDs <- geneIDs[keep,]
data.frame(rownames(counts), rownames(geneIDs))

counts$symbol <- geneIDs$SYMBOL
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]
gene.length <- counts$Length
counts <- counts[,-1]
colnames(counts) <- substr(colnames(counts), start = 13L, stop = 33L)

library(GEOquery)
GSEset <- getGEO("GSE173897", destdir = "./dataset/COVID/HBV_GSE173897/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 46:51)]
rownames(pdata) <- paste(pdata$geo_accession , pdata$title, sep = "_")

keep <- intersect(rownames(pdata), colnames(counts))
pdata <- pdata[keep, ]
counts <- counts[, keep]


library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/HBV_GSE173897/HBV_GSE173897.RData")


###Typhoid Fever cohort ####
counts <- read.delim("./dataset/COVID/TyphoidFever_GSE217667/GSE217667_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE217667", destdir = "./dataset/COVID/TyphoidFever_GSE217667/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 18, 40,43)]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/TyphoidFever_GSE217667/TyphoidFever_GSE217667.RData")


###PBMC night shift nutrition GSE225493####
counts <- read.delim("./dataset/COVID/PBMC_NightShiftNutrition_GSE225493/GSE225493_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE225493", destdir = "./dataset/COVID/PBMC_NightShiftNutrition_GSE225493/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 41,42,44)]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/PBMC_NightShiftNutrition_GSE225493/PBMC_NightShiftNutrition_GSE225493.RData")



###PBMC COVID SA GSE169687####
counts <- read.delim("./dataset/COVID/COVID_SA_GSE169687/GSE169687_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE169687", destdir = "./dataset/COVID/COVID_SA_GSE169687/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 46:53)]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/COVID_SA_GSE169687/PBMC_COVID_SA_GSE169687.RData")



###new-onset T1D GSE124400####
counts <- read.delim("./dataset/COVID/T1D_GSE124400/GSE124400_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE124400", destdir = "./dataset/COVID/T1D_GSE124400/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 50:56)]
pdata <- rbind(pdata, pData(phenoData(GSEset[[2]]))[,c(2,1, 50:56)])

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/T1D_GSE124400/T1D_GSE124400.RData")



###T1D_MS_GSE60424####
counts <- read.delim("./dataset/COVID/T1D_MS_GSE60424/GSE60424_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE60424", destdir = "./dataset/COVID/T1D_MS_GSE60424/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 50:65)]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/T1D_MS_GSE60424/T1D_MS_GSE60424.RData")




###AVSc_GSE218474 whole blood####
counts <- read.delim("./dataset/COVID/AVSc_GSE218474/GSE218474_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE218474", destdir = "./dataset/COVID/AVSc_GSE218474/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 46:49)]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/AVSc_GSE218474/AVSc_GSE218474.RData")


##AMI_GSE249812#####
counts <- read.delim("./dataset/COVID/AMI_GSE249812/GSE249812_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE249812", destdir = "./dataset/COVID/AMI_GSE249812/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 44)]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/AMI_GSE249812/AMI_GSE249812.RData")





##AMI_GSE166780 PBMC#####
counts <- read.delim("./dataset/COVID/AMI_GSE166780/GSE166780_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE166780", destdir = "./dataset/COVID/AMI_GSE166780/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1, 40:42)]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/AMI_GSE166780/AMI_GSE166780.RData")



##STEMI_GSE103182####
counts <- read.delim("./dataset/COVID/STEMI_GSE103182/GSE103182_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

library(GEOquery)
GSEset <- getGEO("GSE103182", destdir = "./dataset/COVID/STEMI_GSE103182/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1)]
pdata$group <- "NSTEMI"
pdata$group[1:15] <- "STEMI"

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]

library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/STEMI_GSE103182/STEMI_GSE103182.RData")





##Sepsis #####
counts <- read.csv("./dataset/COVID/Spesis_PRJCA006118/OMIX006457-02.csv")
genename <- matrix(unlist(strsplit(counts$gene_id, "\\|")), ncol = 2, byrow = T )
genename <- as.data.frame(genename)
counts$gene_id <- genename$V2
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")

rownames(genename) <- genename$V1
rownames(counts) <- genename$V1
anno <- anno[which(!duplicated(anno$EnsemblGeneID)),]
rownames(anno) <- anno$EnsemblGeneID
keep <- intersect(rownames(counts), anno$EnsemblGeneID)
counts <- counts[keep,]
anno <- anno[keep,]
data.frame(rownames(counts), anno$EnsemblGeneID)

counts$gene_id <- anno$Symbol
counts <- aggregate(x = counts[,-1],by = list(counts$gene_id), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]


pdata <- read.csv("./dataset/COVID/Spesis_PRJCA006118/spesis_pdata.csv")
keep <- intersect(pdata$id, colnames(counts))
counts <- counts[,keep]


library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/Spesis_PRJCA006118/Spesis_PRJCA006118.RData")



##Kidney Graft GSE248752######
library(GEOquery)
GSEset <- getGEO("GSE248752", destdir = "./dataset/COVID/KidneyGraft_GSE248752//", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1,38)]

counts <- read.delim("./dataset/COVID/KidneyGraft_GSE248752/GSE248752_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]


library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/KidneyGraft_GSE248752/KidneyGraft_GSE248752.RData")



##Kidney Graft GSE120396######
library(GEOquery)
GSEset <- getGEO("GSE120396", destdir = "./dataset/COVID/KidneyGraft_GSE120396//", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1,11)]

counts <- read.delim("./dataset/COVID/KidneyGraft_GSE120396/GSE120396_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]


library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/KidneyGraft_GSE120396/KidneyGraft_GSE120396.RData")



##Kidney Graft GSE120649######
library(GEOquery)
GSEset <- getGEO("GSE120649", destdir = "./dataset/COVID/KidneyGraft_GSE120649//", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1)]

counts <- read.delim("./dataset/COVID/KidneyGraft_GSE120649/GSE120649_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]


library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/KidneyGraft_GSE120649/KidneyGraft_GSE120649.RData")




##Kidney Graft GSE175718######
library(GEOquery)
GSEset <- getGEO("GSE175718", destdir = "./dataset/COVID/KidneyGraft_GSE175718//", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pdata <- pData(phenoData(GSEset[[1]]))[,c(2,1,42:46)]

counts <- read.delim("./dataset/COVID/KidneyGraft_GSE175718/GSE175718_raw_counts_GRCh38.p13_NCBI.tsv.gz", row.names = 1)
anno <- read.delim("./dataset/COVID/Human.GRCh38.p13.annot.tsv.gz")
data.frame(rownames(counts), anno$GeneID)
counts$symbol <- anno$Symbol
counts <- aggregate(x = counts[,-ncol(counts)],by = list(counts$symbol), FUN = median)
rownames(counts) <- counts$Group.1
counts <- counts[,-1]
counts[1:10,1:10]

keep <- intersect(pdata$geo_accession, colnames(counts))
pdata <- pdata[keep,]


library(edgeR)
y <- DGEList(counts = counts[,], samples = pdata
)
y <- calcNormFactors(y, method = "TMM")
effectiveLibSizes(y )
y$samples$norm.factors
y$samples$norm.factors* y$samples$lib.size
#all.tmm <- sweep(all.count, MARGIN = 2, y$samples$norm.factors, '*')
cpm_by_group_TMM <- cpm(y, normalized.lib.sizes = TRUE)
all.tmm <- round(cpm_by_group_TMM, 6)

load("./dataset/COVID/gene_length.rds")
keep <- intersect(names(gene.length), rownames(all.tmm))
gene.length <- gene.length[keep]
all.tmm <- all.tmm[keep, ]
counts <- counts[keep, ]
#RPKM
rpkm <- edgeR::rpkm(counts, gene.length = gene.length)

#TPM
#TPM = (((mean transcript length in kilobases) x RPKM) / sum(RPKM all genes)) * 10^6.
tmp <- sweep(all.tmm, 1, (gene.length/1000), `/`)
tmp <- sweep(tmp, 2, colSums(tmp), `/`)
tpm <- tmp * 1e6
tpm <- round(tpm, 4)

save(tpm, counts, rpkm, pdata, file = "./dataset/COVID/KidneyGraft_GSE175718/KidneyGraft_GSE175718.RData")





