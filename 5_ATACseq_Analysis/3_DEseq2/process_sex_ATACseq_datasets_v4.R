setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/ATAC-seq_analysis/DEseq2/')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
library('bitops')
library('sva')
library('limma')
library(RColorBrewer)
library(fields)
library(iasva)
library(irlba) 

# process ATAC-dataset
# use Diffbind count on meta peaks

####################################   Neutrophils   ####################################
# read in count matrix
my.ntph1 <- read.csv('../Diffbind/2020-04-18_ATAC_Sex_Aging_Neutrophils_Normalized_count_matrix.txt', sep = "\t", header = T)
my.ntph1$PeakName <- paste(my.ntph1$seqnames,my.ntph1$start,my.ntph1$end,sep = ":")

my.ntph <- my.ntph1[,c(26,6:25)]
colnames(my.ntph) <- c("PeakName" ,
                       "Female_4m_101" , "Female_4m_102" ,"Female_4m_103" ,"Female_4m_104" , "Female_4m_105" ,
                       "Male_4m_111" ,"Male_4m_112" ,"Male_4m_113" ,"Male_4m_114" ,"Male_4m_115" ,
                       "Female_20m_106" ,"Female_20m_107" ,"Female_20m_108" ,"Female_20m_109" ,"Female_20m_110" ,
                       "Male_20m_116" ,"Male_20m_117" ,"Male_20m_118" ,"Male_20m_119" ,"Male_20m_120" )

# get HOMER annotations
my.peak.annot <- read.csv('../Diffbind/HOMER_2020-04-18_ATAC_Sex_Aging_Neutrophils_diffbind_peaks.xls', sep = "\t", header = T)
colnames(my.peak.annot)[1] <- "PeakID"
my.peak.annot$PeakName <- paste(my.peak.annot$Chr,my.peak.annot$Start-1,my.peak.annot$End,sep = ":")

# round counts (DESeq needs integers)
my.ntph[,2:21] <- round(my.ntph[,2:21])
rownames(my.ntph) <- my.ntph$PeakName

# get the peaks with no reads out
my.good <- rownames(my.ntph)[apply(my.ntph[,-1] > 0, 1, sum) > 10 ] # see deseq2 vignette, need to remove too low genes
my.filtered.matrix <- my.ntph[my.good,-1] # 95926 peaks

# color-code
my.colors <- rep("deeppink",15)
my.colors[grep("Female_20m",colnames(my.filtered.matrix))] <- "deeppink4"
my.colors[grep("Male_4m",colnames(my.filtered.matrix))]    <- "deepskyblue"
my.colors[grep("Male_20m",colnames(my.filtered.matrix))]   <- "deepskyblue4"


####################################################################################################################################################
# 1. Run SVA to remove unwanted variation
my.Sex        <- c(rep("F",5),rep("M",5),rep("F",5),rep("M",5))
my.Age        <- c(rep(4,10),rep(20,10))
my.FriP       <- c(0.53,0.30,0.74,0.68,0.59,0.31,0.54,0.67,0.69,0.54,0.80,0.62,0.32,0.58,0.22,0.70,0.27,0.63,0.55,0.44)
my.dup.rates  <- c(18.47,40.40,16.50,14.53,13.26,12.66,32.39,12.35, 9.21,11.02,17.09,18.17,12.25,10.32,12.32,13.05,19.82,11.06,12.88,13.94)

# build design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix ),
                         age  = my.Age,
                         sex  = my.Sex,
                         FriP = my.FriP,
                         dup  = my.dup.rates)

# Set null and alternative models (ignore batch)
mod1 = model.matrix(~ age + sex + FriP + dup, data = dataDesign)
n.sv.be = num.sv(my.filtered.matrix,mod1,method="be") # 2

# apply SVAseq algortihm
my.svseq = svaseq(as.matrix(my.filtered.matrix), mod1, n.sv=n.sv.be, constant = 1)

# remove FrIP, dup and SV, preserve age and sex
my.clean <- removeBatchEffect(log2(my.filtered.matrix + 1),
                              batch=NULL,
                              covariates=cbind(my.svseq$sv,my.FriP,my.dup.rates),
                              design=mod1[,1:3])

# delog and round data for DEseq2 processing
my.filtered.sva <- round(2^my.clean)
write.table(my.filtered.sva, file = paste(Sys.Date(),"Neutrophils_ATAC_postSVA_Counts.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)


####################################################################################################################################################
# 2. DESeq2 on cleaned data
my.outprefix <- paste(Sys.Date(),"Neutrophils_ATAC_DESeq2_Analysis_post_SVA",sep="_")

# design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.sva ),
                         age = my.Age,
                         sex = my.Sex)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.sva,
                              colData = dataDesign,
                              design = ~ age + sex)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# plot dispersion
my.disp.out <- paste(my.outprefix,"dispersion_plot.pdf",sep="_")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix,"MDS_plot.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.05,0.05),
     ylim = c(-0.04,0.08),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 1)
dev.off()

# PCA analysis
my.pos.var <- apply(tissue.cts,1,var) > 0
my.pca <- prcomp(t(tissue.cts[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix,"PCA_plot.pdf",sep="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, pch = 16, col = my.colors,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     #xlim = c(-120,120),
     #ylim = c(-100,100),
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()


# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2, outline = F)
dev.off()

# plot X and Y peak depths expression
pdf(paste(my.outprefix,"_Normalized_counts_Xist_ATAC_barplot.pdf"))
barplot(apply(tissue.cts[rownames(tissue.cts) %in% my.peak.annot$PeakName[my.peak.annot$Gene.Name %in% "Xist"],],2,sum), ylab = "Normalized log2(counts)  X chromosome", las = 2, col = my.colors)
dev.off()

pdf(paste(my.outprefix,"_Normalized_counts_Y_chromosome_ATAC_barplot.pdf"))
barplot(apply(tissue.cts[rownames(tissue.cts) %in% my.peak.annot$PeakName[my.peak.annot$Chr %in% "chrY"],],2,sum), ylab = "Normalized log2(counts) Y chromosome", las = 2, col = my.colors)
dev.off()

###############################################################################################
## a. model aging with sex as covariate  %%%%%%%%%%%%%%
res.age <- results(dds.deseq, name= "age")

### get the heatmap of aging changes at FDR5; exclude NA
res.age <- res.age[!is.na(res.age$padj),]

genes.aging <- rownames(res.age)[res.age$padj < 0.05]
my.num.aging <- length(genes.aging) # 9035

my.heatmap.out <- paste(my.outprefix,"AGING_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Aging significant (FDR<5%), ", my.num.aging, " peaks",sep="")
pheatmap(tissue.cts[genes.aging,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.age, file = paste(Sys.Date(),"Neutrophils_Aging_BOTH.RData", sep ="_"))

###############################################################################################
## b. sex with age as covariate
res.sex <- results(dds.deseq, contrast = c("sex","F","M")) # FC in females over Males

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.sex <- res.sex[!is.na(res.sex$padj),]

genes.sex <- rownames(res.sex)[res.sex$padj < 0.05]
my.num.sex <- length(genes.sex) # 738

my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Sex significant (FDR<5%), ", my.num.sex, " peaks",sep="")
pheatmap(tissue.cts[genes.sex,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.sex, file = paste(Sys.Date(),"Neutrophils_SEX.RData", sep ="_"))


# output result tables of combined analysis to text files
my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_DEseq2_SVA.txt", sep = "_")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.age <- paste(my.outprefix,"AGING_all_genes_statistics.txt", sep = "_")
my.out.stats.sex <- paste(my.outprefix,"SEX_DIM_all_genes_statistics.txt", sep = "_")
write.table(res.age, file = my.out.stats.age , sep = "\t" , row.names = T, quote=F)
write.table(res.sex, file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.age <- paste(my.outprefix,"AGING_FDR5_genes_statistics.txt", sep = "_")
my.out.fdr5.sex <- paste(my.outprefix,"SEX_DIM_FDR5_genes_statistics.txt", sep = "_")
write.table(res.age[genes.aging,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote=F)
write.table(res.sex[genes.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)
################################################################################################


################################################################################################
# annotate Peaks
library(bitops)

##### Log2 normalized counts
tissue.cts <- cbind(rownames(tissue.cts),tissue.cts)
colnames(tissue.cts)[1] <- "PeakName"
tissue.cts.annot <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],tissue.cts))

my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_DEseq2_SVA_PeakAnnot.txt", sep = "_")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

##### Output SEX bed
res.sex <- cbind(rownames(res.sex),res.sex)
colnames(res.sex)[1] <- "PeakName"
res.sex.annot <- data.frame(merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],res.sex))

my.out.stats.sex <- paste(my.outprefix,"SEX_DIM_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.annot, file = my.out.stats.sex , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.sex <- paste(my.outprefix,"SEX_DIM_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.sex.annot[res.sex.annot$padj < 0.05,], file = my.out.fdr5.sex, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.sex.BED <- paste0(my.outprefix,"_SEX_DIM_Background.bed")
write.table(res.sex.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.sex.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.BED <- paste0(my.outprefix,"_SEX_DIM_FDR5.bed")
write.table(res.sex.annot[res.sex.annot$padj < 0.05,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.BED.up <- paste0(my.outprefix,"_SEX_DIM_FDR5_UP.bed")
write.table(res.sex.annot[bitAnd(res.sex.annot$padj < 0.05,res.sex.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.sex.BED.dwn <- paste0(my.outprefix,"_SEX_DIM_FDR5_DWN.bed")
write.table(res.sex.annot[bitAnd(res.sex.annot$padj < 0.05,res.sex.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.sex.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)

##### Output AGE bed
res.age <- cbind(rownames(res.age),res.age)
colnames(res.age)[1] <- "PeakName"

res.age.annot <- merge(my.peak.annot[,c("PeakName","Chr", "Start", "End", "Annotation", "Distance.to.TSS", "Gene.Name","Gene.Description" )],res.age)

my.out.stats.age <- paste(my.outprefix,"AGING_all_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.annot, file = my.out.stats.age , sep = "\t" , row.names = F, quote=F)

my.out.fdr5.age <- paste(my.outprefix,"AGING_FDR5_genes_statistics_PeakAnnot.txt",sep = "_")
write.table(res.age.annot[res.age.annot$padj < 0.05,], file = my.out.fdr5.age, sep = "\t" , row.names = F, quote=F)

my.out.bckgd.age.BED <- paste0(my.outprefix,"_AGING_Background.bed")
write.table(res.age.annot[,c("Chr", "Start", "End","PeakName")], file = my.out.bckgd.age.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.BED <- paste0(my.outprefix,"_AGING_FDR5.bed")
write.table(res.age.annot[res.age.annot$padj < 0.05,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.BED , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.BED.up <- paste0(my.outprefix,"_AGING_FDR5_UP_with_Age.bed")
write.table(res.age.annot[bitAnd(res.age.annot$padj < 0.05,res.age.annot$log2FoldChange >0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.BED.up , sep = "\t" , row.names = F, col.names = F, quote=F)

my.out.fdr5.age.BED.dwn <- paste0(my.outprefix,"_AGING_FDR5_DWN_with_Age.bed")
write.table(res.age.annot[bitAnd(res.age.annot$padj < 0.05,res.age.annot$log2FoldChange <0)>0,c("Chr", "Start", "End","PeakName")], file = my.out.fdr5.age.BED.dwn , sep = "\t" , row.names = F, col.names = F, quote=F)
################################################################################################


#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()
