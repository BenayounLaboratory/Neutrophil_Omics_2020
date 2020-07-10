setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/HOMER_repeats/')
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

# run HOMER repeat counts

####################################   Neutrophils   ####################################
# read in count matrix
my.ntph.rep <- read.csv('HOMER_repeats_raw_matrix.txt', sep = "\t", header = T)
my.ntph.rna <- read.csv('HOMER_RNA_raw_matrix.txt', sep = "\t", header = T)

colnames(my.ntph.rep)[1] <- "TranscriptID"
colnames(my.ntph.rna)[1] <- "TranscriptID"

my.ntph1 <- rbind(my.ntph.rep,my.ntph.rna)
my.ntph2 <- my.ntph1[,-c(2:8)]

colnames(my.ntph2) <- c("TranscriptID" ,
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

my.ntph <- my.ntph2[,c("TranscriptID" ,
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

rownames(my.ntph) <- my.ntph$TranscriptID

# get the genes with no reads out
my.good <- which(apply(my.ntph[,-1]>0, 1, sum) >= 6) # see deseq2 vignette, need to remove too low genes
my.filtered.matrix <- my.ntph[my.good,-1]            

my.genes   <- rownames(my.filtered.matrix)[c(grep('NM',rownames(my.filtered.matrix)), grep('NR',rownames(my.filtered.matrix)))]
my.repeats <- rownames(my.filtered.matrix)[-c(grep('NM',rownames(my.filtered.matrix)), grep('NR',rownames(my.filtered.matrix)))]

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
write.table(my.filtered.sva, file = paste(Sys.Date(),"Neutrophils_HOMER_Repeats_postSVA_Counts.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
########################################################


####################################################################################################################################################
# 2. DESeq2 on cleaned data
my.outprefix <- paste(Sys.Date(),"Neutrophils_DESeq2_Analysis_HOMER_Repeats_post_SVA",sep="_")

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

tissue.cts.repeats <- tissue.cts[my.repeats,]
        
# color-code 
my.colors <- rep("deeppink",16)
my.colors[grep("Female_20m",colnames(tissue.cts.repeats))] <- "deeppink4"
my.colors[grep("Male_4m",   colnames(tissue.cts.repeats))] <- "deepskyblue"
my.colors[grep("Male_20m",  colnames(tissue.cts.repeats))] <- "deepskyblue4"

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts.repeats,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix,"Repeats_ONLY_MDS_plot.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.008, 0.008),
     ylim = c(-0.006,0.006),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las=1)
dev.off()

# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts.repeats,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

###############################################################################################
## a. model aging with sex as covariate  %%%%%%%%%%%%%%
res.age <- results(dds.deseq, name= "age")

### get the heatmap of aging changes at FDR5; exclude NA
res.age <- res.age[!is.na(res.age$padj),]

genes.aging  <- rownames(res.age)[res.age$padj < 0.05]
rep.aging    <- intersect(genes.aging,my.repeats)
my.num.aging <- length(rep.aging)

my.heatmap.out <- paste(my.outprefix,"AGING_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Aging significant (FDR<5%), ", my.num.aging, " Repeats",sep="")
pheatmap(tissue.cts[rep.aging,],
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

genes.sex  <- rownames(res.sex)[res.sex$padj < 0.05]
rep.sex    <- intersect(genes.sex,my.repeats)
my.num.sex <- length(rep.sex)

my.heatmap.out <- paste(my.outprefix,"SEX_DIM_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Sex significant (FDR<5%), ", my.num.sex, " Repeats",sep="")
pheatmap(tissue.cts[rep.sex,],
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

my.out.stats.age <- paste(my.outprefix,"AGING_all_Repeats_statistics.txt",sep = "_")
my.out.stats.sex <- paste(my.outprefix,"SEX_DIM_all_Repeats_statistics.txt",sep = "_")
write.table(res.age[intersect(rownames(res.age),my.repeats),], file = my.out.stats.age , sep = "\t" , row.names = T, quote=F)
write.table(res.sex[intersect(rownames(res.sex),my.repeats),], file = my.out.stats.sex , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.age <- paste(my.outprefix,"AGING_FDR5_Repeats_statistics.txt",sep = "_")
my.out.fdr5.sex <- paste(my.outprefix,"SEX_DIM_FDR5_Repeats_statistics.txt",sep = "_")
write.table(res.age[rep.aging,], file = my.out.fdr5.age, sep = "\t" , row.names = T, quote=F)
write.table(res.sex[rep.sex,], file = my.out.fdr5.sex, sep = "\t" , row.names = T, quote=F)
################################################################################################

#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()

