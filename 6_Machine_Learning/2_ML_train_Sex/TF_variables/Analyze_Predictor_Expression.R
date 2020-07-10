setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/ML_run/TF_predictor_analysis')
options(stringsAsFactors = F)

library(pheatmap)
library(bitops)

# get top predictive TFs

# read feature importance table
my.var.imps <- read.table('../Testing_accuracy/2020-05-20_Variable_Importance_Parsing_RF_GBM_CARET.txt', sep = "\t", header = T)

# variables are in the original order, so the non-TFs are at the top of the table
my.var.imps$VarType <- c(rep("ATAC",2)                    ,
                         rep("Genomic",5)                 ,
                         rep("TF", (nrow(my.var.imps)-7) ) )
  
# Rank of Rank Product
my.sort.rp <- sort(my.var.imps$RANK_Product, index.return = T, decreasing = F)
my.var.imps$Rank_of_Rank_Product[my.sort.rp$ix]    <- 1:nrow(my.var.imps)

# get TFs in top 20
my.top.tfs <- my.var.imps[bitAnd(my.var.imps$Rank_of_Rank_Product <= 20, my.var.imps$VarType %in% "TF" )>0, ]

# Read in normalized gene expression
my.rna.exp  <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/2020-04-02_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt', sep = "\t", header = T)

pdf(paste0(Sys.Date(),"_Top_TF_SexDim_predictor_Gene_Expression.pdf"), onefile = F)
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

my.rna.test <- read.table('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/2020-04-02_Neutrophils_DESeq2_Analysis_STAR_post_SVA_SEX_DIM_all_genes_statistics.txt', sep = "\t", header = T)
my.tf.test <- my.rna.test[my.top.tfs$Feature,]


my.tf.order <- c("Irf8", "Foxm1", "Taf7","Six5", "Tcof1")

pdf(paste0(Sys.Date(),"_Top_TF_SexDim_predictor_Gene_Expression_ORDER.pdf"), onefile = F)
pheatmap(my.rna.exp[my.tf.order,],
         cluster_cols = F,
         cluster_rows = F,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T,
         scale="row",
         main ="Top predicting TF",
         cellwidth = 15,
         cellheight = 10)
dev.off()

