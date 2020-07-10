setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/')
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


####################################   Neutrophils   ####################################
# read in count matrix
my.ntph1 <- read.table('2019-02-04_Aging_BM_Neutrophils_cohort_1_counts.txt', sep = "\t", header = T)
my.ntph2 <- my.ntph1[,-c(2:6)]

colnames(my.ntph2) <- c("GeneName" ,
                       "Female_4m_F1a" ,
                       "Female_4m_F1b" ,
                       "Female_4m_F2a" ,
                       "Female_4m_F2b" ,
                       "Male_4m_M1a"   ,
                       "Male_4m_M1b"   ,
                       "Male_4m_M2a"   ,
                       "Male_4m_M2b"   ,
                       "Female_20m_F3a",
                       "Female_20m_F3b",
                       "Female_20m_F4a",
                       "Female_20m_F4b",
                       "Male_20m_M3a"  ,
                       "Male_20m_M3b"  ,
                       "Male_20m_M4a"  ,
                       "Male_20m_M4b")

my.ntph <- my.ntph2[,c("GeneName" ,
                       "Female_4m_F1a" ,
                       "Female_4m_F1b" ,
                       "Female_4m_F2a" ,
                       "Female_4m_F2b" ,
                       "Female_20m_F3a",
                       "Female_20m_F3b",
                       "Female_20m_F4a",
                       "Female_20m_F4b",
                       "Male_4m_M1a"   ,
                       "Male_4m_M1b"   ,
                       "Male_4m_M2a"   ,
                       "Male_4m_M2b"   ,
                       "Male_20m_M3a"  ,
                       "Male_20m_M3b"  ,
                       "Male_20m_M4a"  ,
                       "Male_20m_M4b")]

my.RINs <- c(6.7, 8, 6.7, 7.2, 6, 7.8, 4.9, 5.7, 6, 4.9, 5.9, 6.9, 6.4, 6.4, 6.1, 6.2)
my.Sex  <- c(rep("F",8),rep("M",8))
my.Age  <- c(rep(4,4),rep(20,4),rep(4,4),rep(20,4))

rownames(my.ntph) <- my.ntph$GeneName

# get the genes with no reads out
my.good <- which(apply(my.ntph[,-1]>0, 1, sum) >= 6) # see deseq2 vignette, need to remove too low genes
my.filtered.matrix <- my.ntph[my.good,-1] # 14784  genes

####################################################################################################################################################
# 1. Run SVA to remove unwanted variation
# build design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                         age = my.Age,
                         sex = my.Sex,
                         rin = my.RINs)

# Set null and alternative models (ignore batch)
mod1 = model.matrix(~ age + sex + rin, data = dataDesign)
n.sv.be = num.sv(my.filtered.matrix,mod1,method="be") # 2

# apply SVAseq algortihm
my.svseq = svaseq(as.matrix(my.filtered.matrix), mod1, n.sv=n.sv.be, constant = 0.1)

# remove RIN and SV, preserve age and sex
my.clean <- removeBatchEffect(log2(my.filtered.matrix + 0.1), 
                              batch=NULL, 
                              covariates=cbind(my.svseq$sv,my.RINs),
                              design=mod1[,1:3])

# delog and round data for DEseq2 processing
my.filtered.sva <- round(2^my.clean-0.1)
write.table(my.filtered.sva, file = paste(Sys.Date(),"Neutrophils_STAR_postSVA_Counts.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
########################################################


####################################################################################################################################################
# 2. DESeq2 on cleaned data
my.outprefix <- paste(Sys.Date(),"Neutrophils_DESeq2_Analysis_STAR_post_SVA",sep="_")

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

# color-code 
my.colors <- rep("deeppink",16)
my.colors[grep("Female_20m",colnames(tissue.cts))] <- "deeppink4"
my.colors[grep("Male_4m",colnames(tissue.cts))]    <- "deepskyblue"
my.colors[grep("Male_20m",colnames(tissue.cts))]   <- "deepskyblue4"

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
     xlim = c(-0.04,0.04),
     ylim = c(-0.025,0.025),
     cex.lab = 1.5,
     cex.axis = 1.5)
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
     xlim = c(-120,120),
     ylim = c(-100,100),
     cex.lab = 1.5,
     cex.axis = 1.5) 
dev.off()


# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

# plot Xist and Ddx3y expression
pdf(paste(my.outprefix,"_Normalized_counts_Xist_expression_barplot.pdf"))
barplot(tissue.cts["Xist",], ylab = "Normalized log2(counts) Xist expression", las = 2, col = my.colors)
dev.off()

pdf(paste(my.outprefix,"_Normalized_counts_Ddx3y_expression_barplot.pdf"))
barplot(tissue.cts["Ddx3y",], ylab = "Normalized log2(counts) Ddx3y expression", las = 2, col = my.colors)
dev.off()


################################################################################################
### a. model age in each sex separately
# get output file prefixes
my.outprefix.age_M <- paste(Sys.Date(),"Neutrophils_DESeq2_model_with_AGING_Males",sep="_")
my.outprefix.age_F <- paste(Sys.Date(),"Neutrophils_DESeq2_model_with_AGING_Females",sep="_")

# get matrix using age as a modeling covariate in FEMALES
dds.F <- DESeqDataSetFromMatrix(countData = my.filtered.sva[, my.Sex %in% 'F'],
                                colData = dataDesign[my.Sex %in% 'F',],
                                design = ~ age)

# run DESeq normalizations and export results
dds.deseq.F <- DESeq(dds.F)
res.age.F <- results(dds.deseq.F, name= "age")

# get matrix using age as a modeling covariate in MALES
dds.M <- DESeqDataSetFromMatrix(countData = my.filtered.sva[, my.Sex %in% 'M'],
                                colData = dataDesign[my.Sex %in% 'M',],
                                design = ~ age)

# run DESeq normalizations and export results
dds.deseq.M <- DESeq(dds.M)
res.age.M <- results(dds.deseq.M, name= "age")


colnames(res.age.F) <- paste(colnames(res.age.F),"F",sep = "_")
colnames(res.age.M) <- paste(colnames(res.age.M),"M",sep = "_")
my.merged.FM <- cbind(res.age.F,res.age.M)
my.merged.FM <- my.merged.FM[!is.na(my.merged.FM$padj_F),]
my.merged.FM <- my.merged.FM[!is.na(my.merged.FM$padj_M),]

my.spear.cor <- cor.test(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M, method = 'spearman')
my.rho <- signif(my.spear.cor$estimate,3)
my.pval <- signif(my.spear.cor$p.value,3)

#### commonly regulated genes
my.F_M.5 <- bitAnd(my.merged.FM$padj_F < 0.05, my.merged.FM$padj_M < 0.05) > 0

pdf(paste(Sys.Date(),"Neutrophils_aging_male_vs_female_FC_scatterplot_FDR5.pdf", sep = "_"))
# par(mar = c(5,4,4,5) + .1)
smoothScatter(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M,
              xlab = "log2(FC per m) in Females with aging",
              ylab = "log2(FC per m) in Males with aging",
              colramp = colorRampPalette(c("white", brewer.pal(9,"YlOrBr"))),
              #postPlotHook = fudgeit,
              cex.lab = 1.5,
              cex.axis = 1.5,
              ylim = c(-0.3,0.3),
              xlim = c(-0.3,0.3)
)

abline(0, 1,  col = "red", lty = "dashed")
abline(h = 0, col = "grey", lty = "dashed")
abline(v = 0, col = "grey", lty = "dashed")
text(0.3, -0.25, paste("Rho =", my.rho), pos = 2)
text(0.3, -0.28, paste("pval =", my.pval), pos = 2)
points(my.merged.FM$log2FoldChange_F[my.F_M.5],my.merged.FM$log2FoldChange_M[my.F_M.5], cex = 1, pch = 16, col = "grey70")
dev.off()

write.table(my.merged.FM[my.F_M.5,], file = paste(Sys.Date(),"Neutrophils_AGING_SexAgreement_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

#### divergently regulated genes
my.F_notM <- bitAnd(my.merged.FM$padj_F < 0.05, my.merged.FM$padj_M > 0.15) > 0
my.M_notF <- bitAnd(my.merged.FM$padj_M < 0.05, my.merged.FM$padj_F > 0.15) > 0

pdf(paste(Sys.Date(),"Neutrophils_aging_male_vs_female_FC_scatterplot_divergent_FDR5.pdf", sep = "_"))
smoothScatter(my.merged.FM$log2FoldChange_F,my.merged.FM$log2FoldChange_M,
              xlab = "log2(FC per m) in Females with aging",
              ylab = "log2(FC per m) in Males with aging",
              colramp = colorRampPalette(c("white", brewer.pal(9,"YlOrBr"))),
              cex.lab = 1.5,
              cex.axis = 1.5,
              ylim = c(-0.3,0.3),
              xlim = c(-0.3,0.3)
)
abline(0, 1,  col = "red", lty = "dashed")
abline(h = 0, col = "grey", lty = "dashed")
abline(v = 0, col = "grey", lty = "dashed")
points(my.merged.FM$log2FoldChange_F[my.F_notM],my.merged.FM$log2FoldChange_M[my.F_notM], cex = 0.5, pch = 16, col = "orchid")
points(my.merged.FM$log2FoldChange_F[my.M_notF],my.merged.FM$log2FoldChange_M[my.M_notF], cex = 0.5, pch = 16, col = "royalblue")
legend("topleft",
       c(paste("M not F,", sum(my.M_notF)),paste("F not M,", sum(my.F_notM))),
       col = c("royalblue","orchid"), pch = 16, pt.cex = 0.5,bty = 'n')
dev.off()

write.table(my.merged.FM[my.M_notF,], file = paste(Sys.Date(),"Neutrophils_AGING_Male.NOT.Female_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(my.merged.FM[my.F_notM,], file = paste(Sys.Date(),"Neutrophils_AGING_Female.NOT.Male_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

save(res.age.F, file = paste(Sys.Date(),"Neutrophils_Female_Aging.RData", sep ="_"))
save(res.age.M, file = paste(Sys.Date(),"Neutrophils_Male_Aging.RData", sep ="_"))
################################################################################################


###############################################################################################
## b. model aging with sex as covariate  %%%%%%%%%%%%%%
res.age <- results(dds.deseq, name= "age")

### get the heatmap of aging changes at FDR5; exclude NA
res.age <- res.age[!is.na(res.age$padj),]

genes.aging <- rownames(res.age)[res.age$padj < 0.05]
my.num.aging <- length(genes.aging)

my.heatmap.out <- paste(my.outprefix,"AGING_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Aging significant (FDR<5%), ", my.num.aging, " genes",sep="")
pheatmap(tissue.cts[genes.aging,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.age, file = paste(Sys.Date(),"Neutrophils_Aging_BOTH.RData", sep ="_"))

###############################################################################################
## c. sex with age as covariate
res.sex <- results(dds.deseq, contrast = c("sex","F","M")) # FC in females over Males

### get the heatmap of sex dimorphic changes at FDR5; exclude NA
res.sex <- res.sex[!is.na(res.sex$padj),]

genes.sex <- rownames(res.sex)[res.sex$padj < 0.05]
my.num.sex <- length(genes.sex)

my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Sex significant (FDR<5%), ", my.num.sex, " genes",sep="")
pheatmap(tissue.cts[genes.sex,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.sex, file = paste(Sys.Date(),"Neutrophils_SEX.RData", sep ="_"))


# output result tables of combined analysis to text files
my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_DEseq2_SVA.txt",sep = "_")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.age <- paste(my.outprefix,"AGING_all_genes_statistics.txt",sep = "_")
my.out.stats.sex <- paste(my.outprefix,"SEX_DIM_all_genes_statistics.txt",sep = "_")
write.table(res.age, file = my.out.stats.age , sep = "\t" , row.names = T, quote=F)
write.table(res.sex, file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.age <- paste(my.outprefix,"AGING_FDR5_genes_statistics.txt",sep = "_")
my.out.fdr5.sex <- paste(my.outprefix,"SEX_DIM_FDR5_genes_statistics.txt",sep = "_")
write.table(res.age[genes.aging,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote=F)
write.table(res.sex[genes.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)
################################################################################################

#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()

