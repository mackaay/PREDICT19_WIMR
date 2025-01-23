library(pheatmap)
library(RColorBrewer)
load(file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./dataset/COVID/gene_length.rds")
load(file = "./COVIDcohort.RData")
load(file = "./dataset/COVID/RF_meanGini.rds")
load(file = "./dataset/COVID/RF_meanGini2.rds")
niteration <- 1220

load(file = "./dataset/COVID/ABIS_RF_sigMat_adjPMN.rds")
load(file = "./dataset/COVID/all.adjusted_PMN.RData")
all.tpm <- adjusted.tpm
load(file = "./dataset/COVID/RF_meanGini_adjPMN.rds")
load(file = "./dataset/COVID/RF_meanGini2_adjPMN.rds")
load(file = "./COVIDcohort.RData")


load(file = "./COVIDcohort.RData")
load(file = "./dataset/COVID/all_adjusted.RData")
load(file = "./dataset/COVID/RF_meanGini_adj.rds")
load( file = "./dataset/COVID/SigMatrix_adj.rds")
all.tpm <- adjusted.tpm
rownames(all.tpm) <- gsub("-", "_", rownames(all.tpm))
niteration <- 2280
load(file = "./dataset/COVID/pal.celltype.rds")



###CIBERSORT methods  WIMR cohort####
library(CIBERSORT)
#sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
#results <- cibersort(sig_matrix, mixture_file)
rownames(tpm) <- gsub("-", "_", rownames(tpm))

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(tpm))
tpm.decon <- tpm[keep,]
decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon), QN= F)



decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.cibersort

library(ggpubr)
manifest <- read.csv("./dataset/COVID/LancetMicrobeSuppl.csv")
rownames(manifest) <- manifest$SampleID

covid <- read.csv("./dataset/COVID/PREDICT-19_AdultCOVIDCohort_ClinicalDatav2_TP_110722.csv")
health <- read.csv("./dataset/COVID/PREDICT-19_HealthyControls_120822.csv")
keep <- intersect(colnames(covid), colnames(health))
covid <- covid[,keep]
health <- health[,keep]
wimr <- rbind(health, covid)
rownames(wimr) <- wimr$SampleID

keep <- intersect(rownames(wimr), rownames(manifest))
wimr <- wimr[keep, ]
manifest <- manifest[keep,]
data.frame(rownames(wimr), rownames(manifest))
rownames(wimr) <- manifest$GEO.accession

#'GSM6730932', 'GSM7507685', 'GSM7507850', 'GSM7507851', 'GSM7507859' 
pData <- pData[which(!duplicated(pData$V2)),]
rownames(pData) <- pData$V2
keep <- intersect(rownames(pData), rownames(wimr))
pData <- pData[keep, ]
wimr <- wimr[keep,]
data.frame(rownames(pData), rownames(wimr))
wimr <- cbind(pData, wimr)
wimr$SRA <- pData$V1

keep <- intersect(rownames(decon.cibersort), wimr$SRA)
decon.all.df <- t(decon.cibersort[keep,1:(ncol(decon.cibersort)-3)])
decon.1.df <- decon.all.df[rowMeans(decon.all.df)>=0.01,]
mat_col <- data.frame(group = wimr$Illness...Control..COVID..Sepsis., 
                      Severity = as.character( wimr$WHO_severity) )
mat_col$Severity <- factor(mat_col$Severity, levels = c("0","1","2","3","4","5","6","7","9"))
rownames(mat_col) <- colnames(decon.all.df)
# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2], 
                   Severity = brewer.pal(9, "YlOrRd"))
names(mat_colors$group) <- unique(wimr$Illness...Control..COVID..Sepsis.)
names(mat_colors$Severity) <- c("0","1","2","3","4","5","6","7","9")
pheatmap(decon.all.df[,order(wimr$WHO_severity)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )


#blood sars cov2 mapping 
wimr$sars_cov2_mapping <- "No"
sars_cov2_mapping <- read.csv("./dataset/COVID/sars_cov2_mapping.csv", header = F)
sars_cov2_mapping <- matrix(as.character(unlist(strsplit(sars_cov2_mapping$V1, "_"))), ncol = 9, byrow = TRUE)
sars_cov2_mapping <- sars_cov2_mapping[,1]

sars_cov2_mapping <- intersect(sars_cov2_mapping, wimr$SRA)
wimr$sars_cov2_mapping[match(sars_cov2_mapping, wimr$SRA)] <- "Yes"

mat_col <- data.frame(group = wimr$Illness...Control..COVID..Sepsis., 
                      Severity = as.character( wimr$WHO_severity), 
                      Sar2_mapping = wimr$sars_cov2_mapping)
rownames(mat_col) <- colnames(decon.all.df)
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2], 
                   Severity = brewer.pal(9, "YlOrRd"), 
                   Sar2_mapping = brewer.pal(8, "Accent")[3:4]
)
names(mat_colors$group) <- unique(wimr$Illness...Control..COVID..Sepsis.)
names(mat_colors$Severity) <- c("0","1","2","3","4","5","6","7","9")
names(mat_colors$Sar2_mapping) <- c("Yes", "No")
pheatmap(decon.1.df[,order(wimr$WHO_severity)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )
rowMeans(decon.1.df)

colnames(wimr)
wimr$bacteria
wimr[125:150,grep("bacter", colnames(wimr))]
wimr["GSM6731002",]
idx <- which(rowSums(wimr[,grep("bacteria_site__", colnames(wimr))]) >0 & is.na(wimr[,"bacteria_site_o"]))
wimr[idx,grep("bacteri*", colnames(wimr))]
idx <- which(rowSums(wimr[,grep("bacterial_species___", colnames(wimr))]) >0 & is.na(wimr[,"bacteria_site_o"]))
wimr[idx,grep("bacteri*", colnames(wimr))]
nrow(wimr[idx,grep("bacteri*", colnames(wimr))])
wimr[,"meta_pneu"]
 

wimr$bacteria_2nd <- "No"
wimr$bacteria_2nd[which(is.na(wimr$bacteria))] <- NA
wimr$bacteria_2nd[idx] <- "Yes"

wimr$diabetes_tmp <- wimr$diabetes
wimr$diabetes[which(wimr$diabetes == 1)] <- "Yes"
wimr$diabetes[which(wimr$diabetes == 0)] <- "No"


mat_col <- data.frame(group = wimr$Illness...Control..COVID..Sepsis., 
                      Severity = as.character( wimr$WHO_severity), 
                      Sar2_mapping = wimr$sars_cov2_mapping, 
                      Bacteria = wimr$bacteria_2nd, 
                      Diabetes = wimr$diabetes)
rownames(mat_col) <- colnames(decon.all.df)
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2], 
                   Severity = brewer.pal(9, "YlOrRd"), 
                   Sar2_mapping = brewer.pal(8, "Accent")[3:4], 
                   Bacteria = brewer.pal(8, "Set1")[3:4], 
                   Diabetes = brewer.pal(8, "Dark2")[3:4])
names(mat_colors$Bacteria) <- c("Yes", "No")
names(mat_colors$group) <- unique(wimr$Illness...Control..COVID..Sepsis.)
names(mat_colors$Severity) <- c("0","1","2","3","4","5","6","7","9")
names(mat_colors$Sar2_mapping) <- c("Yes", "No")
names(mat_colors$Diabetes) <- c("No", "Yes")
pheatmap(decon.1.df[,order(wimr$WHO_severity)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort CIBERSORT meanGini 1220" )
pheatmap(decon.1.df[,order(wimr$bacteria_2nd)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort CIBERSORT meanGini2 2000" )




pheatmap(decon.all.df[,], 
         cluster_cols = T,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )
pheatmap(decon.1.df[,], 
         cluster_cols = T,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )


wimr <- cbind(wimr, t(decon.all.df))
save(wimr, file = "./dataset/COVID/WIMR_decon_cohort.rds")

save(wimr, file = "./dataset/COVID/WIMR_cohort.rds")

save(wimr, file = "./dataset/COVID/WIMR_cohort_2025.rds")



###Japan COVID cohort PBMC ####
load(file = "./dataset/COVID/PBMC_JapanCOVID_GSE206263.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))
library(CIBERSORT)
#sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
#results <- cibersort(sig_matrix, mixture_file)

sig.mat_pbmc <- sig.mat[,grep("Basophils|Erythroblast|matureNeutrophil|Megakaryocyte", colnames(sig.mat), invert = T)]
keep <- intersect(rownames(sig.mat_pbmc[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]
decon.cibersort <- cibersort(sig.mat_pbmc[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon))
rowSums(decon.cibersort[,1:(ncol(decon.cibersort)-3)])
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

pdata$group <- factor(pdata$group, levels = c("Healthy", "Exposed", "Mild", "Severe"))
idx <- grep("0|NA",pdata$days)
mat_col <- data.frame(Days = pdata$days[idx], 
                      Group = pdata$group[idx])
rownames(mat_col) <- colnames(decon.df[,idx])
mat_colors <- list(Group = brewer.pal(4, "Set1"), 
                   Days = brewer.pal(7, "YlOrRd")    )
names(mat_colors$Group) <- unique(pdata$group[idx])
names(mat_colors$Days) <- c("NA","0") #, "3","7","14","21","28")
decon.df.idx <- decon.df[,idx]
pheatmap(decon.df.idx[,order(mat_col$Group)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID Japan cohort CIBERSORT meanGini 1220" )








###COVID Neutrophil cohort GSE212041##### 
#
tpm <- read.delim("./dataset/COVID/GSE212041_Neutrophil_RNAseq_TPM_Matrix.txt.gz")
rownames(tpm) <- gsub("-", "_", rownames(tpm))
#tpm <- aggregate(x = tpm[,3:ncol(tpm)],by = list(tpm$Symbol), FUN = median)
#rownames(tpm) <- tpm$Group.1
#counts <- counts[,-1]
#counts[1:10,1:10]

tpm.decon <- tpm[which(tpm$Symbol %in% rownames(sig.mat[rownames(meanGini)[1:niteration],])),]
rownames(tpm.decon) <- tpm.decon$Symbol
tpm.decon <- tpm.decon[,3:ncol(tpm.decon)]
keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration],]), rownames(tpm.decon))
tpm.decon <- tpm.decon[keep,]
decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,1:20]))

decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

pheatmap(decon.df[,], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         #annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID Neutrophil cohort CIBERSORT meanGini 1220" )




###COVID WBC GSE215865 Mount Sinai (no clinical file)#### 
#longitunial data, need pdata info
library(CIBERSORT)
load(file = "./dataset/COVID/WBC_COVID_GSE215865.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]
tpm.decon <- tpm.decon[,which(!colSums(tpm.decon) == 0) ]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]))
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

pheatmap(decon.df[,], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         #annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID Mount Sanai cohort CIBERSORT meanGini 1220" )

decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]
pheatmap(decon.1.df[,], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         #annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID Mount Sanai cohort CIBERSORT meanGini 1220" )



pdata$gender <- pdata$`Sex:ch1`
pdata$days <- as.numeric( pdata$`days since_first_sample:ch1`)
pdata$covid19 <- pdata$`covid-19 positive:ch1`

####confusion matrix within plates ####


###German COVID cohort (empty)####



###COVID SA PBMC GSE169687 long covid####
#https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-021-02228-6
library(CIBERSORT)
load( file = "./dataset/COVID/COVID_SA_GSE169687/PBMC_COVID_SA_GSE169687.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))

sig.mat_pbmc <- sig.mat[,grep("Basophils|Erythroblast|matureNeutrophil|Megakaryocyte", colnames(sig.mat), invert = T)]

keep <- intersect(rownames(sig.mat_pbmc[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat_pbmc[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]))
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

pdata$`disease severity:ch1` <- factor(pdata$`disease severity:ch1`, levels = c("Healthy", "Mild", "Moderate", "Severe", "Critical"))
pdata$`timepoint:ch1` <- factor(pdata$`timepoint:ch1`, levels = c("Control", "12wpi", "16wpi", "24wpi"))
mat_col <- data.frame(Timepoint = pdata$`timepoint:ch1`, 
                      #PID = pdata$`subject_id:ch1`, 
                      Severity = pdata$`disease severity:ch1`)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Severity = brewer.pal(5, "Set1"), 
                   Timepoint = brewer.pal(4, "YlOrRd")  
                   )
names(mat_colors$Timepoint) <- unique(pdata$`timepoint:ch1`)
names(mat_colors$Severity) <-unique(pdata$`disease severity:ch1`)

decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]
pheatmap(decon.1.df[,order(pdata$`disease severity:ch1`, decreasing = F)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "HBV cohort CIBERSORT meanGini 1220" )



colnames(pdata)[c(5,9,10)] <- c("Severity", "PID", "Timepoint")
pdata <- cbind(pdata, decon.cibersort[,1:(ncol(decon.cibersort)-3)])
library(reshape2)
pdata <- melt(pdata[,4:ncol(pdata)])
library(ggplot2)
library(ggpubr)
pdata[grep("Bex|CD4Teff|CD4Tnaive|CD8Tcm|CD8Teff|Mono_C|Mono_NC|PMN|Tfh|Th1$", pdata$variable), ]%>%
  ggplot( aes(x=Timepoint, y=value)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(aes(color=Severity), size=1, alpha=0.9 ) +scale_color_manual(values=c("#85C1E9", "#E7B800", "#FC4E07", "#BC71F3", "darkgreen"))+
  facet_wrap("variable")+
  theme(legend.position="top",plot.title = element_text(size=11) ) + 
  ggtitle("long COVID SA Cohort") + ylim(c(0,0.38))+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Control", "12wpi"), c("Control", "16wpi"), c("Control", "24wpi")
                                                     ), y_position = c(0.26, 0.30, 0.34))





###HBV cohort (not a lot of difference)####
library(CIBERSORT)
load(file = "./dataset/COVID/HBV_GSE173897/HBV_GSE173897.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]))
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

#pdata$`viral load group:ch1` <- factor(pdata$`viral load group:ch1`, levels = c("<2k", "2-20K", ">20K" ))
pdata$`viral load:ch1` <- as.numeric(pdata$`viral load:ch1`)
mat_col <- data.frame(Status = pdata$`hbv status:ch1`, 
                      vLoad = pdata$`viral load group:ch1`, 
                      Ethnicity = pdata$`ethnicity:ch1`)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Status = brewer.pal(3, "Set1")[1:2], 
                   vLoad = brewer.pal(3, "YlOrRd")    )
names(mat_colors$Status) <- unique(pdata$`hbv status:ch1`)
names(mat_colors$vLoad) <-unique(pdata$`viral load group:ch1`)

decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]
pheatmap(decon.1.df[,order(pdata$`viral load:ch1`, decreasing = F)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "HBV cohort CIBERSORT meanGini 1220" )


pdata <- cbind(pdata, decon.cibersort[,1:(ncol(decon.cibersort)-3)])
colnames(pdata)[3:8] <- c("Ethnicity", "Gender" , "Status", "Treatment", "viral_group", "Load")
pdata$log_load <- log10(pdata$Load)
library(ggpubr)
library(ggfortify)
library(dplyr)
pdata %>%
  ggscatter(x = "log_load", y = "PMN",
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
  )+ stat_cor(method = "pearson", label.x = 6, label.y = 0.025, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "Flow (Percentage)") + xlim(0 , 10) + ylim(0, 0.1)

pdata.pca <- pdata[,9:38]
pdata.pca <- pdata.pca[,which(!colSums(pdata.pca) == 0) ]
prcomp_obj <- prcomp(pdata.pca,  scale. = T)  
autoplot(prcomp_obj, data = pdata.pca, #colour = 'DEPTH',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 5)



pdata$viral_group <- factor(pdata$viral_group, levels = c("<2K", "2-20K", ">20K"))
pdata%>%
  ggplot( aes(x=viral_group, y=PMN, fill=viral_group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#FFEDA0", "#FEB24C", "#F03B20"))+
  #facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("2-20K", ">20K"), 
                                                     c("<2K", ">20K")), y_position = c(0.09, 0.1))
pdata$viral_group <- factor(pdata$viral_group, levels = c("<2K", "2-20K", ">20K"))
pdata%>%
  ggplot( aes(x=viral_group, y=Th1, fill=viral_group)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Th1 ") + scale_fill_manual(values=c("#FFEDA0", "#FEB24C", "#F03B20"))+
  #facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("2-20K", ">20K"), 
                                                     c("<2K", ">20K")), y_position = c(0.25, 0.27))


###Typhoid Fever cohort (different time points )####
#https://www.jci.org/articles/view/169676/figure/1
library(CIBERSORT)
load(file = "./dataset/COVID/TyphoidFever_GSE217667/TyphoidFever_GSE217667.RData")
rownames(tpm) <- gsub("-", "_", rownames(tpm))

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]))
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])



pdata$description <- gsub("not_challenged", "notchallenged", pdata$description)
pdata <- cbind(pdata, as.data.frame(matrix(as.character( unlist(strsplit(pdata$description, split = "_")) ), ncol = 3, byrow = T) ))
colnames(pdata)[6:8] <- c("group", "disease", "day")
pdata$day <- factor(pdata$day, levels = c("V0", "V1", "V7", "D0", "D0.12h", "D7", "TD", "D14"))
pdata$disease <- factor(pdata$disease, levels = c("TD", "nTD", "notchallenged"))


mat_col <- data.frame(Treatment = pdata$`treatment:ch1`,
                      Disease = pdata$disease, 
                      Days = pdata$day)
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Treatment = brewer.pal(4, "Set1"),
                   Disease = brewer.pal(3, "Dark2"), 
                   Days = brewer.pal(8, "Set3"))
names(mat_colors$Treatment) <- unique(pdata$`treatment:ch1`)
names(mat_colors$Disease) <- unique(pdata$disease)
names(mat_colors$Days) <- unique(pdata$day)
decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]
decon.1.df <- decon.1.df[,order(pdata$`treatment:ch1`)]
pdata <- pdata[match(colnames(decon.1.df), rownames(pdata)), ]

pheatmap(decon.1.df, 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Typhoid Fever cohort CIBERSORT meanGini 1220" )
unique(pdata$description)

pdata <- cbind(pdata, t(decon.df))
library(reshape2)
library(ggplot2)
pal.celltype <- brewer.pal(8, "Set1")
pal.celltype <- colorRampPalette(pal.celltype)(29) 

for (i in 1:length(unique(pdata$description))) {
  tmp <- pdata[grep(unique(pdata$description)[i], pdata$description), c("geo_accession", "description", rownames(decon.df))]
  tmp <- melt(tmp)
  
  png(paste0(paste0("./plots/14PercentageBar_TyphoidFever_", unique(pdata$description)[i]), ".png"), units = "cm", width = 10, height = 5, res = 300)
  tmp[,] %>%
    ggplot(data = ., mapping = aes(x = geo_accession, y = value, fill = variable)) +
    geom_col() +
    #geom_text(mapping = aes(label = paste0(Percentage*100, "%") ),              # converting the values to percent
    #          size = 4,                                             # size of the font
    #          position = position_stack(vjust = 0.5)) +             # positioning in the middle
    #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
    scale_fill_manual(values = pal.celltype)+
    #scale_fill_manual(values = c("#0FCC39", "#0FCC98", "#0FC5CC", "#0F96CC", "#0F66CC" ,"#0F37CC" ,"#160FCC" ,"#450FCC", "#CC0FA5", "#CC0F46", "#CC360F", "darkgrey", "orange", "darkorange"))+
    #facet_wrap(.~description) +
    labs(x = "ID",                                              # labelling x axis
         y = "Percentage",                                        # labeling y axis
         title = unique(pdata$description)[i],        # title
         fill = "Cell") +                               # legend
    scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
    theme(
      axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                                 vjust = 0.5),                      # adjusting the position
      axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
      axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
      plot.title = element_text(hjust = 0.5),                       # positioning the plot title
      legend.title = element_text(face = "bold")                    # face the legend title
    ) #+ 
  #geom_hline(yintercept=0.43, linetype="dashed",  color = "darkred", size=2)
  dev.off()
}

i = 2





###Sepsis Spesis_PRJCA006118######
load( file = "./dataset/COVID/Spesis_PRJCA006118/Spesis_PRJCA006118.RData")

library(CIBERSORT)
rownames(tpm) <- gsub("-", "_", rownames(tpm))

keep <- intersect(rownames(sig.mat[rownames(meanGini)[1:niteration ],]), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.cibersort <- cibersort(sig.mat[rownames(meanGini)[1:niteration],], as.matrix(tpm.decon[,]))
decon.cibersort[,1:(ncol(decon.cibersort)-3)] <- round(decon.cibersort[,1:(ncol(decon.cibersort)-3)], 4)
decon.df <- t(decon.cibersort[,1:(ncol(decon.cibersort)-3)])

#pdata <- cbind(pdata, as.data.frame(matrix(as.character( unlist(strsplit(pdata$id, split = "_")) ), ncol = 3, byrow = T) ))
pdata$Treatment[which(pdata$Treatment == 1)] <- "Treated"
pdata$Treatment[which(pdata$Treatment == 0)] <- "NotTreated"

mat_col <- data.frame(Treatment = pdata$Treatment
                      )
rownames(mat_col) <- colnames(decon.df)
mat_colors <- list(Treatment = brewer.pal(3, "Set1"))
names(mat_colors$Treatment) <- unique(pdata$Treatment)

decon.1.df <- decon.df[rowMeans(decon.df)>=0.01,]
decon.1.df <- decon.1.df[,order(pdata$Treatment)]
pdata <- pdata[match(colnames(decon.1.df), pdata$id), ]

pheatmap(decon.1.df, 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Sepsis CIBERSORT meanGini 1220" )
table(pdata$SOFA)





###1a. ABIS signature + NNLS####
library(nnls)
library(pheatmap)
library(RColorBrewer)
abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load(file = "./COVIDcohort.RData")

keep <- intersect(rownames(abis_sig), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
  print(i)
  tmp <- tpm.decon[,i]
  print(colnames(tpm.decon)[i])
  tmp.df <- as.matrix(tpm.decon[,i])
  
  ref.mat <- as.matrix(abis_sig[keep,])
  mod1 <- nnls(ref.mat, tmp.df)
  decon <- mod1$x/sum(mod1$x)
  #colnames(tissue.hyper.mat)
  decon.df <- cbind(decon.df, decon)
  print("merge Done")
}
rownames(decon.df) <- colnames(abis_sig)
decon.all.df <- decon.df
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

rownames(pData) <- pData$V1
pData<- cbind(pData, t(decon.all.df))


# Data frame with column annotations.
mat_col <- data.frame(group = pData$V3)
rownames(mat_col) <- colnames(decon.all.df)
# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2])
names(mat_colors$group) <- unique(pData$V3)
pheatmap(decon.all.df[,order(pData$V3)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )

boxplot(pData$B.Naive ~ pData$V3)
boxplot(pData$B.Memory ~ pData$V3)
boxplot(pData$T.CD4.Naive ~ pData$V3)
boxplot(pData$T.CD4.Memory ~ pData$V3)
boxplot(pData$T.CD8.Naive ~ pData$V3)
boxplot(pData$T.CD8.Memory ~ pData$V3)

goi <- c("DAPP1", "CST3", "FGL2", "GCH1", "CIITA", "UPP1",  "RN7SL1")
pheatmap(tpm[goi,order(pData$V3)], scale = "row", 
         cluster_cols = F, cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","white","darkred"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Gene of Interest" )
tpm[goi,1:10]
counts[goi,1:10]
order(pData$V3)



decon.all.df[,] %>%
  ggplot(data = ., mapping = aes(x = group, y = Percentage, fill = Cluster)) +
  geom_col() +
  geom_text(mapping = aes(label = paste0(Percentage*100, "%") ),              # converting the values to percent
            size = 5,                                             # size of the font
            position = position_stack(vjust = 0.5)) +             # positioning in the middle
  #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
  scale_fill_manual(values = pal.rename)+
  #facet_grid(.~sex) +
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Percentage (Normal)",        # title
       fill = "Cluster") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) #+ 
dev.off()



###1b. signature + NNLS####
library(nnls)
library(pheatmap)
library(RColorBrewer)
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./COVIDcohort.RData")

keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
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
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

rownames(pData) <- pData$V1
pData<- cbind(pData, t(decon.all.df))


# Data frame with column annotations.
mat_col <- data.frame(group = pData$V3)
rownames(mat_col) <- colnames(decon.all.df)
# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2])
names(mat_colors$group) <- unique(pData$V3)
pheatmap(decon.all.df[,order(pData$V3)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )


###1c. overlap + NNLS####
library(nnls)
library(pheatmap)
library(RColorBrewer)
abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load("./dataset/COVID/PBMC_GEO_sigMat.rds")
load(file = "./COVIDcohort.RData")

genelist <- intersect(genelist, rownames(abis_sig))
sig.mat <- sig.mat[genelist,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
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
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

rownames(pData) <- pData$V1
pData<- cbind(pData, t(decon.all.df))


# Data frame with column annotations.
mat_col <- data.frame(group = pData$V3)
rownames(mat_col) <- colnames(decon.all.df)
# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2])
names(mat_colors$group) <- unique(pData$V3)
pheatmap(decon.all.df[,order(pData$V3)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )



###1d. RF + NNLS####
library(nnls)
library(pheatmap)
library(RColorBrewer)
#abis_sig <- read.delim("./dataset/COVID/ABIS_sigmatrixRNAseq.txt", sep = "\t")
load( file = "./dataset/COVID/ABIS_RF_sigMat.rds")
load(file = "./COVIDcohort.RData")
#load(file = "./dataset/COVID/RF_sig.RData")



#sig.mat <- sig.mat[RF.genelist200,]
keep <- intersect(rownames(sig.mat), rownames(tpm))
tpm.decon <- tpm[keep,]

decon.df <- c()
for (i in 1:ncol(tpm)) {
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
colnames(decon.all.df)  <- colnames(tpm.decon)
decon.all.df <- round(decon.all.df, 4)

rownames(pData) <- pData$V1
pData<- cbind(pData, t(decon.all.df))


# Data frame with column annotations.
mat_col <- data.frame(group = pData$V3)
rownames(mat_col) <- colnames(decon.all.df)
# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2])
names(mat_colors$group) <- unique(pData$V3)
pheatmap(decon.all.df[,order(pData$V3)], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )







###2. COVID pheno data #####
library(ggpubr)
manifest <- read.csv("./dataset/COVID/LancetMicrobeSuppl.csv")
rownames(manifest) <- manifest$SampleID

covid <- read.csv("./dataset/COVID/PREDICT-19_AdultCOVIDCohort_ClinicalDatav2_TP_110722.csv")
health <- read.csv("./dataset/COVID/PREDICT-19_HealthyControls_120822.csv")
keep <- intersect(colnames(covid), colnames(health))
covid <- covid[,keep]
health <- health[,keep]
wimr <- rbind(health, covid)
rownames(wimr) <- wimr$SampleID

keep <- intersect(rownames(wimr), rownames(manifest))
wimr <- wimr[keep, ]
manifest <- manifest[keep,]
data.frame(rownames(wimr), rownames(manifest))
rownames(wimr) <- manifest$GEO.accession

#'GSM6730932', 'GSM7507685', 'GSM7507850', 'GSM7507851', 'GSM7507859' 
pData <- pData[which(!duplicated(pData$V2)),]
rownames(pData) <- pData$V2
keep <- intersect(rownames(pData), rownames(wimr))
pData <- pData[keep, ]
wimr <- wimr[keep,]
data.frame(rownames(pData), rownames(wimr))
wimr <- cbind(pData[,5:ncol(pData)], wimr)
wimr$SRA <- pData$V1


###3. Neutrophil correlation #####
neutrophil <- wimr[grep("\\%", wimr$neut ),]
data.frame(neutrophil$Neutrophils, neutrophil$neut)
neutrophil$neut <- gsub("\\%","", neutrophil$neut)
neutrophil$neut <- as.numeric(neutrophil$neut )
neutrophil$Neutrophils <- neutrophil$Neutrophils*100
ggscatter(neutrophil, x = "Neutrophils", y = "neut",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 75, label.y = 20, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "Flow (Percentage)") + xlim(10,100) + ylim(10,100)


neutrophil <- wimr[grep("\\%", wimr$neut ),]
data.frame(neutrophil$Neutrophils.LD, neutrophil$neut)
neutrophil$neut <- gsub("\\%","", neutrophil$neut)
neutrophil$neut <- as.numeric(neutrophil$neut )
neutrophil$Neutrophils.LD <- neutrophil$Neutrophils.LD*100
ggscatter(neutrophil, x = "Neutrophils.LD", y = "neut",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 30, label.y = 90, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "Flow (Percentage)") + xlim(20,100) + ylim(20,100)


###4. COVID severity ####
keep <- intersect(colnames(decon.all.df), wimr$SRA)
decon.all.df <- decon.all.df[,keep]
decon.1.df <- decon.all.df[rowMeans(decon.all.df)>=0.01,]
mat_col <- data.frame(group = wimr$Illness...Control..COVID..Sepsis., 
                      Severity = as.character( wimr$WHO_severity) )
mat_col$Severity <- factor(mat_col$Severity, levels = c("0","1","2","3","4","5","6","7","9"))
rownames(mat_col) <- colnames(decon.all.df)
# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2], 
                   Severity = brewer.pal(9, "YlOrRd"))
names(mat_colors$group) <- unique(wimr$Illness...Control..COVID..Sepsis.)
names(mat_colors$Severity) <- c("0","1","2","3","4","5","6","7","9")
pheatmap(decon.1.df[,order(wimr$WHO_severity)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )


#blood sars cov2 mapping 
wimr$sars_cov2_mapping <- "No"
sars_cov2_mapping <- read.csv("./dataset/COVID/sars_cov2_mapping.csv", header = F)
sars_cov2_mapping <- matrix(as.character(unlist(strsplit(sars_cov2_mapping$V1, "_"))), ncol = 9, byrow = TRUE)
sars_cov2_mapping <- sars_cov2_mapping[,1]

sars_cov2_mapping <- intersect(sars_cov2_mapping, wimr$SRA)
wimr$sars_cov2_mapping[match(sars_cov2_mapping, wimr$SRA)] <- "Yes"

mat_col <- data.frame(group = wimr$Illness...Control..COVID..Sepsis., 
                      Severity = as.character( wimr$WHO_severity), 
                      Sar2_mapping = wimr$sars_cov2_mapping)
rownames(mat_col) <- colnames(decon.all.df)
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2], 
                   Severity = brewer.pal(9, "YlOrRd"), 
                   Sar2_mapping = brewer.pal(8, "Accent")[3:4]
                   )
names(mat_colors$group) <- unique(wimr$Illness...Control..COVID..Sepsis.)
names(mat_colors$Severity) <- c("0","1","2","3","4","5","6","7","9")
names(mat_colors$Sar2_mapping) <- c("Yes", "No")
pheatmap(decon.1.df[,order(wimr$WHO_severity)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )


colnames(wimr)
wimr$bacteria
wimr[125:150,grep("bacter", colnames(wimr))]
wimr["GSM6731002",]
wimr$bacteria_2nd <- "No"
wimr$bacteria_2nd[which(!is.na(wimr$bacteria_date))] <- "Yes"

mat_col <- data.frame(group = wimr$Illness...Control..COVID..Sepsis., 
                      Severity = as.character( wimr$WHO_severity), 
                      Sar2_mapping = wimr$sars_cov2_mapping, 
                      Bacteria = wimr$bacteria_2nd)
rownames(mat_col) <- colnames(decon.all.df)
mat_colors <- list(group = brewer.pal(3, "Set1")[1:2], 
                   Severity = brewer.pal(9, "YlOrRd"), 
                   Sar2_mapping = brewer.pal(8, "Accent")[3:4], 
                   Bacteria = brewer.pal(8, "Set1")[3:4])
names(mat_colors$Bacteria) <- c("Yes", "No")
names(mat_colors$group) <- unique(wimr$Illness...Control..COVID..Sepsis.)
names(mat_colors$Severity) <- c("0","1","2","3","4","5","6","7","9")
names(mat_colors$Sar2_mapping) <- c("Yes", "No")
pheatmap(decon.1.df[,order(wimr$WHO_severity)], 
         cluster_cols = F,cluster_rows = T, display_numbers = F, color = colorRampPalette(c( "darkblue","aquamarine","yellow"))(100), 
         border_color      = NA,
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "COVID WIMR cohort" )


ggscatter(wimr, x = "Monocytes.C", y = "WHO_severity",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.1, label.y = 15, p.accuracy = 0.001, r.accuracy = 0.01)



###5. Boxplot #####
library(viridis)
wimr %>%
  ggplot( aes(x=Illness...Control..COVID..Sepsis., y=Neutrophils, fill=Illness...Control..COVID..Sepsis.)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("COVID vs Health") +
  xlab("")
wimr %>%
  ggplot( aes(x=Illness...Control..COVID..Sepsis., y=CD4Tnaive, fill=Illness...Control..COVID..Sepsis.)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("COVID vs Health") +
  xlab("")
wimr %>%
  ggplot( aes(x=Illness...Control..COVID..Sepsis., y=CD8Tcm, fill=Illness...Control..COVID..Sepsis.)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("COVID vs Health") +
  xlab("")
wimr %>%
  ggplot( aes(x=Illness...Control..COVID..Sepsis., y=Th2, fill=Illness...Control..COVID..Sepsis.)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("COVID vs Health") +
  xlab("")
wimr %>%
  ggplot( aes(x=Illness...Control..COVID..Sepsis., y=Bsm, fill=Illness...Control..COVID..Sepsis.)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("COVID vs Health") +
  xlab("")


###public data, infection, PBMC ####
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206263
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220682
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263756
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237562



