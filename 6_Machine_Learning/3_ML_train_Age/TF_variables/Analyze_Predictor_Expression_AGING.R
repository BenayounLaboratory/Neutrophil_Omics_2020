setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Age_Expression/ML_run/TF_predictor_analysis/')
options(stringsAsFactors = F)

library(pheatmap)
library(bitops)

# get top predictive TFs

# read feature importance table
my.var.imps <- read.table('../Testing_accuracy/2020-06-02_Variable_Importance_Parsing_RF_GBM_CARET_AGING.txt', sep = "\t", header = T)

# variables are in the original order, so the non-TFs are at the top of the table
my.var.imps$VarType <- c(rep("ATAC",2)                    ,
                         rep("Genomic",4)                 ,
                         rep("TF", (nrow(my.var.imps)-6) ) )
  
# Rank of Rank Product
my.sort.rp <- sort(my.var.imps$RANK_Product, index.return = T, decreasing = F)
my.var.imps$Rank_of_Rank_Product[my.sort.rp$ix]    <- 1:nrow(my.var.imps)

# get TFs in top 20
my.top.tfs <- my.var.imps[bitAnd(my.var.imps$Rank_of_Rank_Product <= 20, my.var.imps$VarType %in% "TF" )>0, ]

# Read in normalized gene expression
my.rna.exp  <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt', sep = "\t", header = T)

pdf(paste0(Sys.Date(),"_Top_TF_AGING_predictor_Gene_Expression.pdf"), onefile = F)
pheatmap(my.rna.exp[my.top.tfs$Feature,],
        cluster_cols = F,
        cluster_rows = T,
        colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
        show_rownames = T,
        scale="row",
        main ="Top predicting TF",
        cellwidth = 15,
        cellheight = 10)
dev.off()

my.rna.test <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA_AGING_all_genes_statistics.txt', sep = "\t", header = T)
my.tf.test <- my.rna.test[my.top.tfs$Feature,]
# baseMean log2FoldChange       lfcSE       stat       pvalue         padj
# Xrn2    353.00405  -0.0254821279 0.006033617 -4.2233588 2.406882e-05 0.0004651878
# Klf1    106.01578  -0.0031669049 0.010210446 -0.3101632 7.564368e-01 0.8602227178
# Jarid2  458.83286  -0.0040202377 0.002936108 -1.3692406 1.709241e-01 0.3202923689
# Mtf2    377.25209  -0.0175879477 0.004686865 -3.7526038 1.750073e-04 0.0019338289
# Stat5a  112.79431   0.0270257208 0.008209257  3.2921030 9.944120e-04 0.0068749859
# Foxo1    45.60267   0.0350313469 0.014926113  2.3469838 1.892608e-02 0.0628428548
# Gabpa   394.50456  -0.0053169127 0.004314963 -1.2322036 2.178730e-01 0.3830591534
# Cbx8     23.39706   0.0167764427 0.012062478  1.3907957 1.642874e-01 0.3110730436
# Gtf2b   161.18948  -0.0081404918 0.005635191 -1.4445813 1.485755e-01 0.2899668007
# Kat2a    70.85780  -0.0200606463 0.008652089 -2.3185899 2.041728e-02 0.0667410999
# Sirt6    12.84383   0.0282497545 0.015549960  1.8167092 6.926165e-02 0.1666828154
# Six5      6.54291   0.0525968284 0.029602934  1.7767438 7.561041e-02 0.1776681132
# Taf7   1168.19800   0.0005781761 0.004633319  0.1247866 9.006925e-01 0.9473864617
# Nod2    106.34638   0.0189546802 0.006458301  2.9349327 3.336200e-03 0.0169552270

my.tf.sig <- c( "Mtf2","Nod2","Stat5a","Xrn2","Foxo1", "Kat2a")

pdf(paste0(Sys.Date(),"_Top_TF_AGING_predictor_Gene_Expression_sig_ORDER.pdf"), onefile = F)
pheatmap(my.rna.exp[my.tf.sig,],
        cluster_cols = F,
        cluster_rows = F,
        colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
        show_rownames = T,
        scale="row",
        main ="Top predicting TF",
        cellwidth = 15,
        cellheight = 10)
dev.off()

