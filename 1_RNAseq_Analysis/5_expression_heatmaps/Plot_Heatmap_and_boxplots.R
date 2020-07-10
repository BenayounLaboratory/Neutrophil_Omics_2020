setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/')
library('grDevices')
library('bitops')
library('pheatmap')

# 2020-04-02
# get granule component heatmap

# 2020-06-18
# replot, and add boxplots for Padi4 and atgl/Pnpla2 

#### load data
load('2020-05-21_Neutrophils_SEX.RData')
res.sex.df <- data.frame(res.sex)
res.sex.df <- cbind(rownames(res.sex.df),res.sex.df)
colnames(res.sex.df)[1] <- "GeneNames"

tissue.cts <- read.table('2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt')


# load gene lists
my.1.gran.v2     <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/Neutrophil_Gene_Lists/Primary_Granules_v2.txt', sep = "\t")
my.2.gran.v2     <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/Neutrophil_Gene_Lists/Secondary_Granules_v2.txt', sep = "\t")
my.3.gran.v2     <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/Neutrophil_Gene_Lists/Tertiary_Granules_v2.txt', sep = "\t")
my.netosis       <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/Neutrophil_Gene_Lists/NETosis_regulation.txt', sep = "\t")

my.var.null <- apply(tissue.cts,1,var)

my.1.gran.v2     <- intersect(unique(my.1.gran.v2$V1) , rownames(tissue.cts)[my.var.null != 0])
my.2.gran.v2     <- intersect(unique(my.2.gran.v2$V1) , rownames(tissue.cts)[my.var.null != 0])
my.3.gran.v2     <- intersect(unique(my.3.gran.v2$V1) , rownames(tissue.cts)[my.var.null != 0])
my.netosis       <- intersect(unique(my.netosis$V1)   , rownames(tissue.cts)[my.var.null != 0])

pdf(paste0(Sys.Date(),"_Primary_Granule_LONG_Gene_Expression.pdf"), onefile = F)
pheatmap(tissue.cts[my.1.gran.v2,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Primary Granules",
         cellwidth = 15,
         cellheight = 10)
dev.off()

pdf(paste0(Sys.Date(),"_Secondary_Granule_LONG_Gene_Expression.pdf"), onefile = F)
pheatmap(tissue.cts[my.2.gran.v2,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Secondary Granules",
         cellwidth = 15,
         cellheight = 10)
dev.off()

pdf(paste0(Sys.Date(),"_Tertiary_Granule_LONG_Gene_Expression.pdf"), onefile = F, height = 10)
pheatmap(tissue.cts[my.3.gran.v2,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Tertiary Granules",
         cellwidth = 15,
         cellheight = 10)
dev.off()

############################################################
pdf(paste0(Sys.Date(),"_Sex-linked_Gene_Expression.pdf"), onefile = F)
pheatmap(tissue.cts[c("Xist","Ddx3y", "Eif2s3y","Uty","Kdm5c","Kdm6a"),],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Sex-linked genes",
         cellwidth = 15,
         cellheight = 10)
dev.off()

############################################################
pdf(paste0(Sys.Date(),"_NETOsis_Gene_Expression.pdf"), onefile = F)
pheatmap(tissue.cts[my.netosis,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="NETosis",
         cellwidth = 15,
         cellheight = 10)
dev.off()

############################################################
pdf(paste0(Sys.Date(),"_Padi4_Gene_Expression.pdf"), onefile = F)
pheatmap(tissue.cts["Padi4",],
         cluster_cols = F,
         cluster_rows = F,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Padi4",
         cellwidth = 15,
         cellheight = 10)
dev.off()

###################################################################################################
######## get boxplots

my.padi4.exp <- data.frame("YF" = as.numeric(tissue.cts["Padi4",1:4]  ) ,
                           "OF" = as.numeric(tissue.cts["Padi4",5:8]  ) ,
                           "YM" = as.numeric(tissue.cts["Padi4",9:12] ) ,
                           "OM" = as.numeric(tissue.cts["Padi4",13:16]) )

pdf(paste0(Sys.Date(),"_Padi4_Gene_Expression_boxplot.pdf"), onefile = F, height = 5, width = 3)
boxplot(my.padi4.exp, 
        col = c("deeppink", "deeppink4","deepskyblue", "deepskyblue4"),
        ylab = "DESeq2 normalized log2 counts",
        ylim = c(10,14)
        )
beeswarm::beeswarm(my.padi4.exp, add = T, pch = 16)
dev.off()



my.Pnpla2.exp <- data.frame("YF" = as.numeric(tissue.cts["Pnpla2",1:4]  ) ,
                           "OF" = as.numeric(tissue.cts["Pnpla2",5:8]  ) ,
                           "YM" = as.numeric(tissue.cts["Pnpla2",9:12] ) ,
                           "OM" = as.numeric(tissue.cts["Pnpla2",13:16]) )

pdf(paste0(Sys.Date(),"_Pnpla2_Atgl_Gene_Expression_boxplot.pdf"), onefile = F, height = 5, width = 3)
boxplot(my.Pnpla2.exp, 
        col = c("deeppink", "deeppink4","deepskyblue", "deepskyblue4"),
        ylab = "DESeq2 normalized log2 counts",
       ylim = c(7,10)
)
beeswarm::beeswarm(my.Pnpla2.exp, add = T, pch = 16)
dev.off()
