library(openxlsx)
library(pheatmap)
library(ggplot2)
library(openxlsx)
setwd("e://Metabolomics2/181010_Mouse_neutrophils_aging_BB/")

data_pHILIC <- read.csv("BB_Mouse_neutrophils_HILIC_pos.csv")
data_nHILIC <- read.csv("BB_Mouse_neutrophils_HILIC_neg.csv")
data_pRPLC <- read.csv("BB_Mouse_neutrophils_RPLC_pos.csv")
data_nRPLC <- read.csv("BB_Mouse_neutrophils_RPLC_neg.csv")

data_pHILIC$Compounds_ID <- paste("pHILIC", data_pHILIC$Compounds_ID, sep = "_")
data_nHILIC$Compounds_ID <- paste("nHILIC", data_nHILIC$Compounds_ID, sep = "_")
data_pRPLC$Compounds_ID <- paste("pRPLC", data_pRPLC$Compounds_ID, sep = "_")
data_nRPLC$Compounds_ID <- paste("nRPLC", data_nRPLC$Compounds_ID, sep = "_")

# Combine datasets
data <- rbind(data_pHILIC,data_nHILIC,data_pRPLC,data_nRPLC)
data <- data[,-1]
data$Mode <- c(rep("pHILIC", nrow(data_pHILIC)), rep("nHILIC", nrow(data_nHILIC)), rep("pRPLC", nrow(data_pRPLC)), rep("nRPLC", nrow(data_nRPLC))) 

rm(list = c("data_pHILIC_QC","data_nHILIC_QC","data_pRPLC_QC","data_nRPLC_QC", "data_pHILIC","data_nHILIC","data_pRPLC","data_nRPLC",
            "df", "tmp"))

##########################################################
#### Data curation ####
###########################################################
data_ID <- data[which(data$Metabolite != 0),]
data_ID <- data_ID[,-grep("QC", colnames(data_ID))]
write.csv(data_ID, "BB_Mouse_Neutrophils_metabolomics.csv")

data_ID <- read.csv("BB_Mouse_Neutrophils_metabolomics.csv", row.names = 1)
colnames(data_ID)

###############################################################################################################################################
############# Hierarchical Clustering #############
###############################################################################################################################################
# Load file
data_nz <- data_ID[,11:26]

# Log transformation
data_nz.log10 <- log10(data_nz)

# Range scaling
data_nz.mean <- apply(data_nz.log10,1,mean)
data_nz.log10.centered <- data_nz.log10 - data_nz.mean
data_nz.max <- apply(data_nz.log10.centered,1,max)
data_nz.min <- apply(data_nz.log10.centered,1,min)
data_nz.range <- data_nz.max - data_nz.min
data_nz.norm <- data_nz.log10.centered/data_nz.range
data_nz.norm <- data_nz.norm[complete.cases(data_nz.norm),]

# Clustering using hcluster from amap package, use spearman and complete 
library('amap')
sample = t(data_nz.norm)

## perform clustering
hc <- hcluster(t(data_nz.norm), method = "spearman", diag = FALSE, upper = FALSE,
               link = "complete", members = NULL, nbproc = 2,
               doubleprecision = TRUE)
plot(hc, hang = -1)

# Calculate PCs
data.pca <- prcomp(t(data_nz.norm),center=F, scale=F)
summary(data.pca)
sum.pca<-summary(data.pca)
imp.pca<-sum.pca$importance;
std.pca<-imp.pca[1,] # standard deviation
var.pca<-imp.pca[2,] # variance explained by each PC
cum.pca<-imp.pca[3,] # cummulated variance explained

# PCA plot 2D
library(ggplot2)
Comp <- c(2,4) # Choose PC to plot
ggdata <- data.frame(data.pca$x[,Comp[1]],data.pca$x[,Comp[2]])
ggdata$Category <- 0
ggdata$Category[grep("F_4m", rownames(ggdata))] <- "Young_Female"
ggdata$Category[grep("M_4m", rownames(ggdata))] <- "Young_Male"
ggdata$Category[grep("F_20m", rownames(ggdata))] <- "Old_Female"
ggdata$Category[grep("M_20m", rownames(ggdata))] <- "Old_Male"
colnames(ggdata) <- c(paste("PC", Comp, sep=""),"Category")
# stat_ellipse is not part of the base ggplot package
source("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R") 
ggplot(ggdata, aes(x=PC2,y=PC4,fill=factor(Category), label=rownames(ggdata))) +
  stat_ellipse(geom="polygon", level=0.95, alpha=0.25, border=NA) +
  geom_point(size=5, shape=21, color="black") +
  guides(color=guide_legend("Category"),fill=guide_legend("Category")) +
  scale_x_continuous(paste("PC2 (", round(100*var.pca[Comp[1]],1),"%)",sep="")) +
  scale_y_continuous(paste("PC7 (", round(100*var.pca[Comp[2]],1),"%)",sep="")) +
  theme_bw() + theme(legend.justification=c(1,1), legend.position=c(1,1)) +
  theme(axis.title = element_text(face="bold", size=16)) + theme(legend.title = element_text(colour="white"))

# Statistical Analysis
stats <- data_ID

for (i in 1:nrow(stats)) {
  y <- as.numeric(stats[i,grep("4m_", colnames(stats))])
  x <- as.numeric(stats[i,grep("20m_", colnames(stats))])
  stats$logFC[i] <- log2((mean(x)/mean(y)))
  tt <- t.test(x,y, paired= F, alternative="two.sided")
  stats$pval[i] <- tt$p.value
}
stats$FDR <- p.adjust(stats$pval, method = "BH")

hist(stats$pval, breaks = 100, col="darkblue")
hist(stats$FDR, breaks = 100, col="darkblue")

length(which(stats$FDR < 0.15))

pheatmap(stats[which(stats$FDR < 0.15),1:16], cluster_rows= T, cluster_cols=T, clustering_distance_cols = "euclidean", scale="row",
         clustering_method = "complete", color=colorRampPalette(c("blue", "white", "red"))(100)) 
