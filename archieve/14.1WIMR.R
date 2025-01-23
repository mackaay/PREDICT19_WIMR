library(fmsb)
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

load(file = "./dataset/COVID/WIMR_cohort.rds")
load(file = "./dataset/COVID/WIMR_cohort.rds")
load(file = "./dataset/COVID/WIMR_cohort_2025.rds")
wimr

##0. gender and age ####
wimr$Severity_Group <- "Healthy"
idx.group <- grep("1|2|3", wimr$WHO_severity)
wimr$Severity_Group[idx.group] <- "Mild" 
idx.group <- grep("4|5|6", wimr$WHO_severity)
wimr$Severity_Group[idx.group] <- "Moderate" 
idx.group <- grep("7|9", wimr$WHO_severity)
wimr$Severity_Group[idx.group] <- "Severe" 

table(wimr$sex, wimr$Severity_Group)
table(wimr$sex, wimr$outcome)

library(ggplot2)
library(ggpubr)
library("scales")

pt <- table(wimr$sex, wimr$Severity_Group)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var2 <- factor(pt$Var2, levels = c("Healthy", "Mild", "Moderate", "Severe"))
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_classic(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = c( "#CD534CFF","#86868699"))+
  xlab("Severity") +
  ylab("Proportion") 



wimr.df <- wimr[, c( "CD8Tcm", "CD8Teff","CD4Teff", "Tfh", "Th1", "Treg","Neutrophil", "PMNMDSC",  "MMDSC",  "Plasmablasts", "Severity_Group",  "sex")]
p1 <- wimr.df[complete.cases(wimr.df),]%>%
    ggplot( aes(x=sex, y=PMNMDSC, fill=sex)) +
    geom_boxplot(outlier.size=-1) +
    geom_jitter(color="black", size=1, alpha=0.5) +
    theme(legend.position="none",plot.title = element_text(size=11) ) + 
    ggtitle("PMN-MDSC ") + scale_fill_manual(values=c( "#CD534CFF","#86868699"))+
    facet_wrap("Severity_Group")+
    ylab("Percentage")+ geom_signif(comparisons = list(c("female", "male")), y_position = 0.4)  
p4 <-  wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=sex, y=CD8Teff, fill=sex)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("CD8Teff ") + scale_fill_manual(values=c( "#CD534CFF","#86868699"))+
  facet_wrap("Severity_Group")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("female", "male")), y_position = 0.18) 
p5 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=sex, y=Tfh, fill=sex)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Tfh ") + scale_fill_manual(values=c( "#CD534CFF","#86868699"))+
  facet_wrap("Severity_Group")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("female", "male")), y_position = 0.25) 
wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=sex, y=MMDSC, fill=sex)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("M-MDSC ") + scale_fill_manual(values=c( "#CD534CFF","#86868699"))+
  facet_wrap("Severity_Group")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("female", "male")), y_position = 0.25) 
p6 <-wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=sex, y=Plasmablasts, fill=sex)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Plasmablasts ") + scale_fill_manual(values=c( "#CD534CFF","#86868699"))+
  facet_wrap("Severity_Group")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("female", "male")), y_position = 0.25) 
p3 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=sex, y=Th1, fill=sex)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Th1 ") + scale_fill_manual(values=c( "#CD534CFF","#86868699"))+
  facet_wrap("Severity_Group")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("female", "male")), y_position = 0.25) 
p2 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=sex, y=Neutrophil, fill=sex)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Neutrophil ") + scale_fill_manual(values=c( "#CD534CFF","#86868699"))+
  facet_wrap("Severity_Group")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("female", "male")), y_position = 0.45) 
ggarrange(p1, p2, p3, p4,p5, p6,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2)  

  

#Age 
ggscatter(wimr, x = "WHO_severity", y = "age",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 7.5, label.y = 90, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Correlation of age and severity",x ="WHO Severity", y = "Age") + xlim(0,10) + ylim(25,100)

ggscatter(wimr, x = "PMN", y = "age",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.4, label.y = 40, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Correlation of age and PMN",x ="PMN proportion", y = "Age") + xlim(0,0.65) + ylim(25,100)
ggscatter(wimr, x = "Th1", y = "age",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.3, label.y = 40, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Correlation of age and Th1",x ="Th1 proportion", y = "Age") + xlim(0,0.45) + ylim(25,100)
ggscatter(wimr, x = "Neutrophil", y = "age",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.5, label.y = 40, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Correlation of age and Neutrophil",x ="mature Neutrophil proportion", y = "Age") + xlim(0,0.55) + ylim(25,100)



##1. Radarchart Severity ####
idx.cell <- c( "CD8Tcm", "CD8Teff","CD4Teff", "Tfh", "Th1", "Treg","Neutrophil", "PMNMDSC",  "Plasmablasts" ,  "MMDSC")
wimr.df <- wimr[,idx.cell]
idx.group <- grep("healthy", wimr$V3)
wimr.chart <- colMeans(wimr.df[idx.group,])
idx.group <- grep("1|2|3", wimr$WHO_severity)
wimr.chart <- rbind(wimr.chart, colMeans(wimr.df[idx.group,]))
idx.group <- grep("4|5|6", wimr$WHO_severity)
wimr.chart <- rbind(wimr.chart, colMeans(wimr.df[idx.group,]))
idx.group <- grep("7|9", wimr$WHO_severity)
wimr.chart <- rbind(wimr.chart, colMeans(wimr.df[idx.group,]))
wimr.chart <- as.data.frame(wimr.chart)
rownames(wimr.chart) <- c("Health", "Mild", "Moderate", "Severe")
max_min <- data.frame(
  CD8Tcm = c(0.25, 0), CD8Teff = c(0.25, 0), CD4Teff = c(0.25, 0), Tfh = c(0.25, 0),Th1 = c(0.25, 0),
  Treg = c(0.25, 0),  Neutrophil = c(0.25, 0), PMNMDSC = c(0.25, 0),
  MMDSC = c(0.25, 0), Plasmablasts = c(0.25, 0)
)
rownames(max_min) <- c("Max", "Min")
wimr.chart <- rbind(max_min, wimr.chart)

radarchart(wimr.chart)
create_beautiful_radarchart(wimr.chart[,], 
                            caxislabels = c(0, 0.5, 0.1, 0.15, 0.2, 0.25),
                            color = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))

colors <- c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3")
titles <- c("Health", "Mild", "Moderate", "Severe")
# Reduce plot margin using par()
# Split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(2,2))
# Create the radar chart
for(i in 1:4){
  create_beautiful_radarchart(
    data = wimr.chart[c(1, 2, i+2), ], caxislabels = c(0, 0.5, 0.1, 0.15, 0.2, 0.25),
    color = colors[i], title = titles[i]
  )
}
par(op)


par(mfrow = c(1,1))
dev.off()





##2. WBC correlation ####
library(ggplot2)
library(ggpubr)

neutrophil <- wimr[grep("\\%", wimr$neut ),]
neutrophil <- neutrophil[,c("neut", "Neutrophil", "PMNMDSC")]
neutrophil$neut <- gsub("\\%","", neutrophil$neut)
neutrophil$neut <- as.numeric(neutrophil$neut )
neutrophil$Neutrophil <- neutrophil$Neutrophil*100
neutrophil$PMN <- neutrophil$PMN*100

tmp <- wimr[grep("[1-9]", wimr$neut ),]
tmp <- tmp[grep("\\%", tmp$neut, invert = T),]
tmp$neut <- as.numeric(tmp$neut )
tmp$wbc <- as.numeric(tmp$wbc )
tmp$neut <- tmp$neut / tmp$wbc
tmp <- tmp[,c("neut", "Neutrophil", "PMNMDSC")]
tmp <- tmp[tmp$neut< 1 ,]
tmp$Neutrophil <- tmp$Neutrophil*100
tmp$PMN <- tmp$PMN*100
tmp$neut <- tmp$neut*100

neutrophil <- rbind(neutrophil, tmp)

neutrophil$Neutrophils <- neutrophil$Neutrophil + neutrophil$PMN
ggscatter(neutrophil, x = "Neutrophils", y = "neut",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 40, label.y = 40, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "Flow (Percentage)") + xlim(0,60) + ylim(0,100)

#WBC
wimr$wbc
wbc <- wimr[grep("[1-9]", wimr$wbc ),"wbc"]
wbc <- cbind(wbc, wimr[grep("[1-9]", wimr$wbc ),c(343:371)]  )
wbc <- wbc[,c(1:12, 14:16, 18:ncol(wbc))]  
wbc$WBC <- rowSums(wbc[,2:ncol(wbc)])

ggscatter(wbc, x = "WBC", y = "wbc",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 40, label.y = 40, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="WBC",x ="RNAseq (Percentage)", y = "Flow (Count x10^9)") + xlim(5,70) + ylim(25,100)


#
wimr$lymp  
lymphocyte <- wimr[grep("\\%", wimr$lymp ),]
lymphocyte <- lymphocyte[,c("lymp", "Bex", "Bnaive", "Bnsm",    "Bsm", "CD4Teff", 
                            "CD4Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem", "CD8Tnaive", 
                            "MAIT", "NK",  "Tfh" , "Th1",   "Th17", "Th2",   "Treg", "Tgd")]
lymphocyte$lymp <- gsub("\\%","", lymphocyte$lymp)
lymphocyte$lymp <- as.numeric(lymphocyte$lymp )
lymphocyte[,2:ncol(lymphocyte)]<- lymphocyte[,2:ncol(lymphocyte)] *100

tmp <- wimr[grep("[1-9]", wimr$lymp ),]
tmp <- tmp[grep("\\%", tmp$lymp, invert = T),]
tmp$lymp <- as.numeric(tmp$lymp )
tmp$wbc <- as.numeric(tmp$wbc )
tmp$lymp <- tmp$lymp / tmp$wbc
tmp <- tmp[,c("lymp", "Bex", "Bnaive", "Bnsm",    "Bsm", "CD4Teff", 
              "CD4Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem", "CD8Tnaive", 
              "MAIT", "NK",  "Tfh" , "Th1",   "Th17", "Th2",   "Treg", "VD2")]
tmp <- tmp[tmp$lymp< 1 ,]
tmp<- tmp*100

lymphocyte <- rbind(lymphocyte, tmp)

lymphocyte$lymphocytes <- rowSums(lymphocyte[,2:ncol(lymphocyte)])
ggscatter(lymphocyte, x = "lymphocytes", y = "lymp",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 25, label.y = 50, p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(title="lymphocyte",x ="RNAseq (Percentage)", y = "Flow (Percentage)") + xlim(20,80) + ylim(0,60)

colMax(lymphocyte)
colMin(lymphocyte)



nlr <- wimr[grep("[1-9]", wimr$lymp ),]
nlr <- nlr[grep("\\%", nlr$lymp, invert = T),]
nlr$lymp <- as.numeric(nlr$lymp )
nlr$neut <- as.numeric(nlr$neut )
nlr$nlr <- nlr$neut / nlr$lymp
nlr <- nlr[,c("nlr","outcome", "Bex", "Bnaive", "Bnsm",    "Bsm", "CD4Teff", 
              "CD4Tnaive", "CD8Tcm", "CD8Teff", "CD8Tem", "CD8Tnaive", 
              "MAIT", "NK",  "Tfh" , "Th1",   "Th17", "Th2",   "Treg", "VD2")]
table(nlr$nlr, nlr$outcome)




##3. Outcomes ####
manifest <- read.csv("./dataset/COVID/LancetMicrobeSuppl.csv")
manifest <- manifest[which(manifest$GEO.accession %in% rownames(wimr)),]
keep <- intersect(manifest$GEO.accession, rownames(wimr))
rownames(manifest) <- manifest$GEO.accession
manifest <- manifest[match(rownames(wimr), manifest$GEO.accession),]
data.frame(manifest$GEO.accession, rownames(wimr))
wimr$V2 <- manifest$SampleID
colnames(wimr)[1:3] <- c("Severity", "SampleID", "Group" )
wimr$Severity <- NA
idx.group <- grep("1|2|3", wimr$WHO_severity)
wimr$Severity[idx.group] <- "Mild"
idx.group <- grep("4|5|6", wimr$WHO_severity)
wimr$Severity[idx.group] <- "Moderate"
idx.group <- grep("7|9", wimr$WHO_severity)
wimr$Severity[idx.group] <- "Severe"

table(wimr$Severity, wimr$outcome)
wimr$Outcome_hospital <- NA
idx.group <- grep("0", wimr$outcome)
wimr$Outcome_hospital[idx.group] <- "Dead"
idx.group <- grep("1", wimr$outcome)
wimr$Outcome_hospital[idx.group] <- "Alive"

wimr.df <- wimr[grep("Dead|Alive", wimr$Outcome_hospital),c("Outcome_hospital", "PMNMDSC", "Severity", 
                                                            "CD8Teff", "Tfh", "Th1", "Neutrophil","CD8Tcm", "cMono","MMDSC",
                                                            "NK", "Treg", "CD4Teff", "Plasmablasts", "pDC")]
wimr.df$Neut_ratio <- wimr.df$PMN/wimr.df$Neutrophil
p1 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=PMNMDSC, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN MDSC") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.4)
p4 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=CD8Teff, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("CD8Teff ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity") + ylim(0, 0.25)+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.2)
p5 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=Tfh, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Tfh ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+ylim(0, 0.35)+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.3)
wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=MMDSC, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("MMDSC ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.4)
wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=NK, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("NK ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.15)
wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=cMono, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("cMono ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.3)
p6 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=Plasmablasts, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Plasmablasts ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.2)
wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=pDC, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("pDC ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.0075)
wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=pDC, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Neut_ratio ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.0075)
p3 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=Th1, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Th1 ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.35)
p2 <- wimr.df[complete.cases(wimr.df),]%>%
  ggplot( aes(x=Outcome_hospital, y=Neutrophil, fill=Outcome_hospital)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Neutrophil ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+ ylim(0, 0.6)+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Alive", "Dead")), y_position = 0.5)
ggarrange(p1, p2, p3, p4,p5, p6,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2)





##4. Longitudinal ####
library(limma)
library(dplyr)
longitudial <- wimr[grep("d[0-9]", wimr$SampleID),]
longitudial$SampleID
tmp <- as.data.frame(strsplit2(longitudial$SampleID, "d"))
colnames(tmp) <- c("tmp_V1", "tmp_V2")
longitudial <- cbind(longitudial, tmp )
n_occur <- data.frame(table(longitudial$tmp_V1)) 
longitudial <- longitudial[which(longitudial$tmp_V1 %in% n_occur$Var1[n_occur$Freq > 1]),]

longitudial$Outcome_hospital <- NA
idx.group <- grep("0", longitudial$outcome)
longitudial$Outcome_hospital[idx.group] <- "Dead"
idx.group <- grep("1", longitudial$outcome)
longitudial$Outcome_hospital[idx.group] <- "Alive"
longitudial$Group <- "Health"
idx.group <- grep("1|2|3", longitudial$WHO_severity)
longitudial$Group[idx.group] <- "Mild"
idx.group <- grep("4|5|6", longitudial$WHO_severity)
longitudial$Group[idx.group] <- "Moderate"
idx.group <- grep("7|9", longitudial$WHO_severity)
longitudial$Group[idx.group] <- "Severe"


longitudial <- longitudial[,c("tmp_V1", "tmp_V2", "PMN", "Neutrophil", "Th1", "CD8Teff", "Tfh", "Plasmablasts","Mono_C","Outcome_hospital", "bacteria_2nd" ,"Group")]
longitudial$tmp_V2 <- as.numeric(longitudial$tmp_V2)
longitudial <- longitudial[order(longitudial$tmp_V1),]
longitudial$Group[4] <- "Severe"
longitudial$Group[10] <- "Severe"
longitudial$Group[22] <- "Severe"
longitudial$Group[76:79] <- "Moderate"
longitudial$Group[91] <- "Moderate"
longitudial$Group[63:66] <- "Moderate"

p1 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=PMN, group=tmp_V1, colour=Outcome_hospital ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p2 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Neutrophil, group=tmp_V1, colour=Outcome_hospital ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p3 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Th1, group=tmp_V1, colour=Outcome_hospital ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p4 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=CD8Teff, group=tmp_V1, colour=Outcome_hospital ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p5 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Tfh, group=tmp_V1, colour=Outcome_hospital ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p6 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Plasmablasts, group=tmp_V1, colour=Outcome_hospital ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Mono_C, group=tmp_V1, colour=Outcome_hospital ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
ggarrange(p1, p2, p3, p4,p5,p6,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 1, nrow = 6)


p1 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=PMN, group=tmp_V1, colour=bacteria_2nd ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p2 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Neutrophil, group=tmp_V1, colour=bacteria_2nd ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p3 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Th1, group=tmp_V1, colour=bacteria_2nd ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p4 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=CD8Teff, group=tmp_V1, colour=bacteria_2nd ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p5 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Tfh, group=tmp_V1, colour=bacteria_2nd ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
p6 <- ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Plasmablasts, group=tmp_V1, colour=bacteria_2nd ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
ggplot(longitudial[complete.cases(longitudial),], aes(x=tmp_V2, y=Mono_C, group=tmp_V1, colour=bacteria_2nd ) ) + 
  geom_line(size=1.5)+ scale_color_manual(values=c("#85C1E9", "#F7DC6F","black"))+ 
  facet_wrap(Group~.)
ggarrange(p1, p2, p3, p4,p5,p6,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 1, nrow = 6)



###5. ICU###
wimr$los_hosp <- as.numeric(wimr$los_hosp)
wimr$los_icu <- as.numeric(wimr$los_icu)
wimr$total_symp 
wimr$total_co 

ggscatter(wimr, x = "PMN", y = "los_hosp", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.1, label.y = 125, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="PMN",x ="RNAseq (Percentage)", y = "Days in hospital") + xlim(0,0.35) + ylim(0, 150)
ggscatter(wimr, x = "Neutrophil", y = "los_hosp", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.1, label.y = 125, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "Days in hospital") + xlim(0,0.35) + ylim(0,250)
ggscatter(wimr, x = "CD8Teff", y = "los_hosp",color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.05, label.y = 125, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="CD8Teff",x ="RNAseq (Percentage)", y = "Days in hospital") + xlim(0,0.2) + ylim(0,150)
ggscatter(wimr, x = "CD8Tcm", y = "los_hosp",color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.05, label.y = 125, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="CD8Tcm",x ="RNAseq (Percentage)", y = "Days in hospital") + xlim(0,0.35) + ylim(0,250)
ggscatter(wimr, x = "Tfh", y = "los_hosp",color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.05, label.y = 125, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="Tfh",x ="RNAseq (Percentage)", y = "Days in hospital") + xlim(0,0.25) + ylim(0,250)
ggscatter(wimr, x = "Plasmablasts", y = "los_hosp",color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.05, label.y = 125, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="Plasmablasts",x ="RNAseq (Percentage)", y = "Days in hospital") + xlim(0,0.25) + ylim(0,250)


ggscatter(wimr, x = "PMN", y = "los_icu",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 50, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="PMN",x ="RNAseq (Percentage)", y = "Days in ICU") + xlim(0,0.35) + ylim(0,250)
ggscatter(wimr, x = "Neutrophil", y = "los_icu",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 50, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "Days in ICU") + xlim(0,0.35) + ylim(0,250)
ggscatter(wimr, x = "CD8Teff", y = "los_icu",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 50, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="CD8Teff",x ="RNAseq (Percentage)", y = "Days in ICU") + xlim(0,0.35) + ylim(0,250)
ggscatter(wimr, x = "CD8Tcm", y = "los_icu",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 50, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="CD8Tcm",x ="RNAseq (Percentage)", y = "Days in ICU") + xlim(0,0.35) + ylim(0,250)



ggscatter(wimr, x = "PMN", y = "total_symp", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  #scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="PMN",x ="RNAseq (Percentage)", y = "Total Symtoms") + xlim(0,0.35) + ylim(0, 15)
ggscatter(wimr, x = "Neutrophil", y = "total_symp", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "No. of Symptoms") + xlim(0,0.35) + ylim(0,15)
ggscatter(wimr, x = "CD8Teff",y = "total_symp", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="CD8Teff",x ="RNAseq (Percentage)", y = "Days in ICU") + xlim(0,0.35) + ylim(0,15)
ggscatter(wimr, x = "CD8Tcm", y = "total_symp", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="CD8Tcm",x ="RNAseq (Percentage)", y = "Days in ICU") + xlim(0,0.35) + ylim(0,15)


ggscatter(wimr, x = "PMN", y = "total_co", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  #scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="PMN",x ="RNAseq (Percentage)", y = "No. of Complications") + xlim(0,0.35) + ylim(0, 15)
ggscatter(wimr, x = "Neutrophil", y = "total_co", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.1, label.y = 12, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+ 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#BC71F3"))+ theme(legend.position = "none")+
  labs(title="Neutrophil",x ="RNAseq (Percentage)", y = "No. of Complications") + xlim(0,0.35) + ylim(0,15)
ggscatter(wimr, x = "CD8Teff",y = "total_co", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="CD8Teff",x ="RNAseq (Percentage)", y = "No. of Complications") + xlim(0,0.35) + ylim(0,15)
ggscatter(wimr, x = "Th1", y = "total_co", color = "Severity_Group",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
)+ stat_cor(method = "pearson", label.x = 0.2, label.y = 10, p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_wrap(Severity_Group ~.)+
  labs(title="Th1",x ="RNAseq (Percentage)", y = "No. of Complications") + xlim(0,0.35) + ylim(0,15)



##5. co-infection #####
wimr$bacteria_2nd
wimr$sars_cov2_pcr_date
wimr$influenza_a[is.na(wimr$influenza_a)] <- 99
wimr$influenza_b[is.na(wimr$influenza_b)] <- 99
wimr$parainfluenza[is.na(wimr$parainfluenza)] <- 99
wimr$rsv[is.na(wimr$rsv)] <- 99
wimr$adenovirus[is.na(wimr$adenovirus)] <- 99
wimr$enterovirus[is.na(wimr$enterovirus)] <- 99
wimr$rhinovirus[is.na(wimr$rhinovirus)] <- 99
wimr$meta_pneu[is.na(wimr$meta_pneu)] <- 99
wimr[,c("influenza_a", "influenza_b", "parainfluenza", "rsv", "adenovirus", "enterovirus", "rhinovirus", "meta_pneu")]

wimr$virus_coinfection <-  wimr$influenza_a + wimr$influenza_b + wimr$parainfluenza + wimr$rsv + wimr$adenovirus + wimr$enterovirus + wimr$rhinovirus + wimr$meta_pneu
wimr$virus_coinfection[which(wimr$virus_coinfection == 792)] <- NA
wimr$virus_coinfection[which(wimr$virus_coinfection > 200)] <- "No"
wimr$virus_coinfection[which(wimr$virus_coinfection == "0" )] <- "No"
wimr$virus_coinfection[which(wimr$virus_coinfection == "1" )] <- "Yes"


wimr.virus <- wimr[which(!is.na(wimr$virus_coinfection)), ]
wimr.virus$Severity_Group
table(wimr.virus$Severity_Group, wimr.virus$virus_coinfection)

wimr.virus[,c("PMN", "virus_coinfection", "Severity")]%>%
  ggplot( aes(x=virus_coinfection, y=PMN, fill=virus_coinfection)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Percentage")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.4)



wimr$bacteria
wimr$bacteria <- gsub(99, NA, wimr$bacteria)
wimr.bac <- wimr[which(!is.na(wimr$bacteria)),]
wimr.bac$bacteria <- gsub("0", "No", wimr.bac$bacteria)
wimr.bac$bacteria <- gsub("1", "Yes", wimr.bac$bacteria)

p1 <- wimr.bac[,]%>%
  ggplot( aes(x=bacteria, y=PMN, fill=bacteria)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.4)
p2 <- wimr.bac[,]%>%
  ggplot( aes(x=bacteria, y=Mono_C, fill=bacteria)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_C ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.25)
p3 <- wimr.bac[,]%>%
  ggplot( aes(x=bacteria, y=Mono_I, fill=bacteria)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_I ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.075)
p4 <- wimr.bac[,]%>%
  ggplot( aes(x=bacteria, y=Mono_NC, fill=bacteria)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_NC ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.15)
p5 <- wimr.bac[,]%>%
  ggplot( aes(x=bacteria, y=Th1, fill=bacteria)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Th1 ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.21)
ggarrange(p1, p2, p3, p4,p5, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 3, nrow = 2)



##6. Vaccination #####
wimr$vaccination

wimr.vaccination <- wimr[which(!is.na(wimr$vaccination)),]
wimr.vaccination$vaccination
wimr.vaccination <- wimr.vaccination[,-2]

wimr.vaccination$vaccination <- ifelse(wimr.vaccination$vaccination == 0 , "No", "Yes")
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=PMN, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.4)

wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Mono_C, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_C ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.25)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Mono_I, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_I ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.075)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Mono_NC, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Mono_NC ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.15)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Th1, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Th1 ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.2)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Bex, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Bex ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.15)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Bnaive, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Bnaive ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.15)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Bnsm, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Bnsm ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.05)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Bsm, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Bsm ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.05)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=CD8Tcm, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("CD8Tcm ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.2)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=CD8Teff, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("CD8Teff ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.2)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=CD8Tem, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("CD8Tem ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.05)
wimr.vaccination%>%
  ggplot( aes(x=vaccination, y=Plasmablasts, fill=vaccination)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("Plasmablasts ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  facet_wrap("Severity_Group")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.2)
#Bex Bnaive Bnsm Bsm CD4Teff CD4Tnaive CD8Tcm CD8Teff CD8Tem






##7. ECMO#####
wimr[which(wimr$ecmo == 1 ), "PMN"]
table(wimr$ecmo, wimr$Severity_Group)
wimr.severe <- wimr[which(wimr$Severity_Group == "Severe"),]
wimr.severe$ecmo <- ifelse(wimr.severe$ecmo==0, "No", "Yes")
wimr.severe <- wimr.severe[,which(!duplicated(colnames(wimr.severe)))]
wimr.severe%>%
  ggplot( aes(x=ecmo, y=PMN, fill=ecmo)) +
  geom_boxplot(outlier.size=-1) +
  geom_jitter(color="black", size=1, alpha=0.5) +
  theme(legend.position="none",plot.title = element_text(size=11) ) + 
  ggtitle("PMN ") + scale_fill_manual(values=c("#85C1E9", "#F7DC6F"))+
  #facet_wrap("Severity")+
  ylab("Proportion")+ geom_signif(comparisons = list(c("Yes", "No")), y_position = 0.21)






