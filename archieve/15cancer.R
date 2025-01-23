library(nnls)
library(pheatmap)
library(RColorBrewer)

###Melanoma PD1 ####
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/gene_length.rds")
load(file = "./dataset/COVID/RF_meanGini.rds")
load(file = "./dataset/COVID/RF_meanGini2.rds")
niteration <- 1220

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

####RF 500, CIBERSORT####
library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, mixture_file)

idx.tissue <- grep("Plasma|Eryth|Mega", colnames(sig.mat), invert = T) # in tissue, removing Plasma cells, Erythro, and Mega cell types
decon.cibersort <- cibersort(sig.mat[rownames(meanGini2)[1:niteration],idx.tissue], as.matrix(tpm))
decon.cibersort <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])
decon.cibersort_mat <- cibersort(sig_matrix, as.matrix(tpm))
#ABIS dataset is suitalbe for PBMC not for tumor
decon.cibersort_mat <- t(decon.cibersort_mat[,1:22])

keep <- intersect(colnames(decon.cibersort), rownames(pData.cancer))
decon.cibersort <- decon.cibersort[,keep]
pData.cancer <- pData.cancer[keep,]
decon.cibersort_mat <- decon.cibersort_mat[,keep]

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
pheatmap(decon.cibersort_mat[,order(pData.cancer$Best.RECIST.response)], 
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

p1 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=PMN, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("PMN PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.15)
p2 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),]  %>%
  ggplot( aes(x=Response, y=CD8Tem, fill=Response)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("CD8Tem PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.075)
p3 <- pData.cancer[grep("PRE", pData.cancer$RNA.Sequencing),]  %>%
  ggplot( aes(x=Response, y=CD8Teff, fill=Response)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("CD8Teff PRE") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.075)
ggarrange(p1, p2, p3, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

p1 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),] %>%
  ggplot( aes(x=Response, y=PMN, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("PMN EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.01)
p2 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),]  %>%
  ggplot( aes(x=Response, y=CD8Tem, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("CD8Tem EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.125)
p3 <- pData.cancer[grep("EDT", pData.cancer$RNA.Sequencing),]  %>%
  ggplot( aes(x=Response, y=CD8Teff, fill=Response)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + 
  ggtitle("CD8Teff EDT") +
  ylab("Percentage")+ geom_signif(comparisons = list(c("noResponse", "Response")), y_position = 0.125)
ggarrange(p1, p2, p3, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

####RF500, NNLS####
keep <- intersect(rownames(sig.mat[,]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm.decon)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(sig.mat)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm)
decon.all.df <- round(decon.all.df, 4)

keep <- intersect(rownames(pData.cancer), colnames(decon.all.df))
pData.cancer <- pData.cancer[keep,]
decon.all.df <- decon.all.df[,keep]
pData.cancer<- cbind(pData.cancer, t(decon.all.df))


mat_col <- data.frame(Response = pData.cancer$Best.RECIST.response, 
                      Time = pData.cancer$RNA.Sequencing)
rownames(mat_col) <- colnames(decon.all.df)
# List with colors for each annotation.
mat_colors <- list(Best.RECIST.response = brewer.pal(4, "Set1"),
                   Time = brewer.pal(3, "Dark2")[1:2])
names(mat_colors$Best.RECIST.response) <- unique(pData.cancer$Best.RECIST.response)
names(mat_colors$Time) <- unique(pData.cancer$RNA.Sequencing)

pheatmap(decon.all.df[,order(pData.cancer$Best.RECIST.response)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Melanoma IPIPD1 cohort" )



###irAE melanoma GSE186143####
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
library(e1071)
tpm <- read.delim("./dataset/COVID/GSE186143_bulk_TPM_samples.txt.gz", row.names = 1)
colnames(tpm) <- gsub("\\.","-", colnames(tpm))
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186143

library(GEOquery)
GSEset <- getGEO("GSE186143", destdir = "./dataset/COVID/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1, 48:53)]
rownames(sampleinfo) <- sampleinfo$`patient:ch1`
data.frame(colnames(tpm), rownames(sampleinfo))


####CIBERSORT####
library(CIBERSORT)
decon.irAE <- cibersort(sig.mat[rownames(meanGini)[1:500],], as.matrix(tpm))
decon.irAE <- round(decon.irAE[,1:(ncol(decon.irAE)-3)] , 4)
decon.irAE <- t(decon.irAE)

#sampleinfo <- cbind(sampleinfo, decon.irAE)


####NNLS####
library(nnls)
keep <- intersect(rownames(sig.mat[,]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.irAE <- c()
for (i in 1:ncol(tpm.decon)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(sig.mat[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.irAE <- cbind(decon.irAE, decon)
  print("merge Done")
}
rownames(decon.irAE) <- colnames(sig.mat)
colnames(decon.irAE)  <- colnames(tpm)
decon.irAE <- round(decon.irAE, 4)


keep <- intersect(rownames(sampleinfo), colnames(decon.irAE))
sampleinfo <- sampleinfo[keep,]
decon.irAE <- decon.irAE[,keep]
#sampleinfo<- cbind(sampleinfo, t(decon.irAE))


library(RColorBrewer)
library(pheatmap)
mat_col <- data.frame(group = sampleinfo$`therapy type:ch1`, 
                      irAE = as.character( sampleinfo$`highest irae grade:ch1`) )
mat_col$irAE <- factor(mat_col$irAE, levels = c("0","1","2","3","4"))
rownames(mat_col) <- colnames(decon.irAE)
mat_colors <- list(group = brewer.pal(3, "Set1")[1:3], 
                   irAE = brewer.pal(5, "YlOrRd"))
names(mat_colors$group) <- unique(sampleinfo$`therapy type:ch1`)
names(mat_colors$irAE) <- c("0","1","2","3","4")
pheatmap(decon.irAE[,order(sampleinfo$`highest irae grade:ch1`)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Melanoma irAE cohort" )

pheatmap(decon.irAE[rowMeans(decon.irAE)>=0.01,order(sampleinfo$`highest irae grade:ch1`)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Melanoma irAE cohort" )








