setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Age_Expression/Features/Feature_Engineering')
options(stringsAsFactors = F)
library(bitops)

source('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/Features/Feature_Engineering/Parsing_gmt_features_functions.R')

# Prep features for aging ML, based on already run Sex ML

###################################################################################################
####### A. Gene Expression feature (response curve)
###################################################################################################

my.expressed.genes <- read.csv('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA_AGING_all_genes_statistics.txt', sep = "\t", header = T)

my.RNA.features              <- data.frame(matrix(0,length(rownames(my.expressed.genes)),3))
rownames(my.RNA.features)    <- rownames(my.expressed.genes)
colnames(my.RNA.features)    <- c("FDR", "logFC", "AGING")

my.RNA.features$FDR   <- my.expressed.genes$padj
my.RNA.features$logFC <- my.expressed.genes$log2FoldChange
my.RNA.features$AGING <- "NOT_AGING"

for (i in 1:nrow(my.RNA.features)) {
  
  if (my.RNA.features$FDR[i] < 0.05) {
    
    if (my.RNA.features$logFC[i] > 0) {
      
      my.RNA.features$AGING[i] <- "UP"
      
    }  else if (my.RNA.features$logFC[i] < 0) {
      
      my.RNA.features$AGING[i] <- "DOWN"
      
    }
  }
}

my.RNA.features  <- data.frame(my.RNA.features[my.RNA.features$AGING != "NOT_AGING",])
my.RNA.features$AGING <- factor(my.RNA.features$AGING)

summary(my.RNA.features$AGING)
# DOWN   UP 
# 1972 1681 

save(my.RNA.features, 
     file = paste(Sys.Date(),"RNA_features_for_ML_AGING.RData", sep = "_"))

###################################################################################################
####### B. TF target features
###################################################################################################

# Harmonizome
my.coll.Chea            <- read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/CHEA_Harmonizome_2020-03-25.gmt")
my.coll.hm.JASPAR       <- read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/JASPAR_TF_targets_Harmonizome_2020-03-31.gmt")
my.coll.hm.ENCODE       <- read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/ENCODE_TFs_Harmonizome_2020-03-25.gmt")

# GEO (Harmonizome and Genome research)
my.coll.GEO             <- read.gmt("../../../ML_Predict_Sex_Dim_Expression/Features/TF_GMT_cleanup/2020-05-20_GEO_Harmonizome_Target_summary.gmt")

#################
# convert to mouse gene syntax rules
my.coll.Chea$gene             <- firstup(tolower( my.coll.Chea$gene          ))           
my.coll.hm.JASPAR$gene        <- firstup(tolower( my.coll.hm.JASPAR$gene     ))     
my.coll.hm.ENCODE$gene        <- firstup(tolower( my.coll.hm.ENCODE$gene     ))     
my.coll.GEO$gene              <- firstup(tolower( my.coll.GEO$gene           ))     

my.coll.Chea$ont             <- firstup(tolower( my.coll.Chea$ont          ))          
my.coll.hm.JASPAR$ont        <- firstup(tolower( my.coll.hm.JASPAR$ont     ))     
my.coll.hm.ENCODE$ont        <- firstup(tolower( my.coll.hm.ENCODE$ont     ))     
my.coll.GEO$ont              <- firstup(tolower( my.coll.GEO$ont           ))     

my.colls.summary             <- rbind(my.coll.Chea      ,
                                      my.coll.hm.JASPAR ,
                                      my.coll.hm.ENCODE ,
                                      my.coll.GEO        )
dim(my.colls.summary)
#2,338,722       2

my.colls.summary <- unique(my.colls.summary)
dim(my.colls.summary)
# 2,141,774       2

# get unique target sets
length(unique(my.colls.summary$ont        ))  # 391

################################
# Parse into a feature matrix
my.expressed.genes <- read.csv('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA_AGING_all_genes_statistics.txt', sep = "\t", header = T)
my.sig.genes <- rownames(my.expressed.genes)[my.expressed.genes$padj < 0.05]
  
# filter on expressed/significant genes
my.summary.sig       <-  my.colls.summary$gene      %in%   my.sig.genes
my.colls.summary.sig <-  my.colls.summary[my.summary.sig,]

# parse into matrix
my.TF.features             <- data.frame(matrix(0,length(my.sig.genes),length(unique(my.colls.summary.sig$ont))))
rownames(my.TF.features)   <- my.sig.genes
colnames(my.TF.features)   <- unique(my.colls.summary.sig$ont)

for (i in 1:nrow(my.colls.summary.sig)) {
  my.TF.features[my.colls.summary.sig$gene[i],my.colls.summary.sig$ont[i]] <- 1
}

dim(my.TF.features) # 3653  391

### Retain TFs with more than 25 targets among significant genes
my.TF.features.ALL <- my.TF.features[,apply(my.TF.features,2,sum) > 25]
dim(my.TF.features.ALL) # 3653  357

# TF genes with detectable expression
my.expressed.TFs <- colnames(my.TF.features.ALL) %in% rownames(my.expressed.genes) # 249
my.TF.features.FINAL <- my.TF.features.ALL[,my.expressed.TFs]

save(my.TF.features.FINAL, file = paste(Sys.Date(),"summary_TF_features_from_GMT_for_ML_AGING.RData", sep = "_"))

library(pheatmap)
pdf(paste0(Sys.Date(),"_correlation_heatmap_TF_target_features_collapsed_AGING.pdf"))
pheatmap(cor(my.TF.features.ALL), show_rownames = F, show_colnames = F)
dev.off()

###################################################################################################
####### C. Sequence features
###################################################################################################

####################################
# Homer promoter CpG/GC content
my.homer.prom.data.1 <- read.csv('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/Features//Constant_features/2020-04-28_HOMER_Neutrophil_Promoter_features.txt', sep = "\t", header = T)
colnames(my.homer.prom.data.1)[1] <- "PeakID"

my.homer.prom.data.2 <- my.homer.prom.data.1[,c("Gene.Name", "Chr", "GC.", "CpG.")]
colnames(my.homer.prom.data.2) <- c("Gene.Name","Chrosomome" ,"Promoter_GC", "Promoter_CpG")

# clean unplaced scaffolds
unique(my.homer.prom.data.2$Chrosomome)
# [1] "chr16"                "chr5"                 "chr17"                "chr2"                
# [5] "chr10"                "chr13"                "chr15"                "chr9"                
# [9] "chr4"                 "chr3"                 "chr11"                "chr7"                
# [13] "chr1"                 "chr6"                 "chr18"                "chr8"                
# [17] "chrX"                 "chr12"                "chr14"                "chr19"               
# [21] "chrY"                 "chr1_GL456221_random" "chr4_GL456216_random" "chrUn_JH584304"      
# [25] "chr1_GL456211_random" "chrX_GL456233_random"

my.homer.prom.data.2$Chrosomome[my.homer.prom.data.2$Chrosomome %in% "chr1_GL456221_random"] <- "chr1"
my.homer.prom.data.2$Chrosomome[my.homer.prom.data.2$Chrosomome %in% "chr4_GL456216_random"] <- "chr4"
my.homer.prom.data.2$Chrosomome[my.homer.prom.data.2$Chrosomome %in% "chrUn_JH584304"]       <- "chrUn"
my.homer.prom.data.2$Chrosomome[my.homer.prom.data.2$Chrosomome %in% "chr1_GL456211_random"] <- "chr1"
my.homer.prom.data.2$Chrosomome[my.homer.prom.data.2$Chrosomome %in% "chrX_GL456233_random"] <- "chrX"

####################################
# Homer gene features
my.homer.gene.data.1 <- read.csv('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/Features/Constant_features/2020-04-28_HOMER_Neutrophil_Gene_features.txt', sep = "\t", header = T)
colnames(my.homer.gene.data.1)[1] <- "PeakID"

my.homer.gene.data.2 <- my.homer.gene.data.1[,c("Gene.Name", "GC.", "CpG.")]
colnames(my.homer.gene.data.2) <- c("Gene.Name", "Exonic_GC", "Exonic_CpG")

#########################
my.expressed.genes <- read.csv('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/DEseq2/STAR/2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA_AGING_all_genes_statistics.txt', sep = "\t", header = T)
my.sig.genes <- rownames(my.expressed.genes)[my.expressed.genes$padj < 0.05]

my.homer.prom.data.sig <- my.homer.prom.data.2[my.homer.prom.data.2$Gene.Name %in% my.sig.genes, ]
my.homer.gene.data.sig <- my.homer.gene.data.2[my.homer.gene.data.2$Gene.Name %in% my.sig.genes, ]

my.seq.feat.alls   <- merge(my.homer.prom.data.sig, my.homer.gene.data.sig, by = "Gene.Name")

save(my.seq.feat.alls, 
     file = paste(Sys.Date(),"HOMER_sequence_features_for_ML_AGING.RData", sep = "_"))

dim(my.seq.feat.alls)
###[1] 3421    7
# only 3421 genes could have features extracted

###################################################################################################
####### D. ATAC-seq features
###################################################################################################

my.ATAC.seq.data <- read.csv('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/Features/ATAC/HOMER_Neutrophil_ATAC_seq_reads.txt', sep = "\t", header = T)
colnames(my.ATAC.seq.data)[1] <- "PeakID"

### rename to short
my.ATAC.seq.data.2 <- my.ATAC.seq.data[,c(16,20:23)]
colnames(my.ATAC.seq.data.2) <- c("Gene.Name"                       ,                                                 
                                  "Neutrophil_ATACseq_4m_F_Merged"  ,
                                  "Neutrophil_ATACseq_4m_M_Merged"  ,
                                  "Neutrophil_ATACseq_21m_F_Merged" ,
                                  "Neutrophil_ATACseq_21m_M_Merged"  )
rownames(my.ATAC.seq.data.2) <- my.ATAC.seq.data.2$Gene.Name

my.ATAC.seq.data.3  <- data.frame("Gene.Name"                      =  my.ATAC.seq.data.2$Gene.Name,
                                  "ATAC_TSS_FPKM_average"          =  apply(my.ATAC.seq.data.2[,-1],1,mean),
                                  "log2_ATAC_TSS_FPKM_old_vs_yg"   =  log2( (apply(my.ATAC.seq.data.2[,c(2,3)],1,mean)/ (apply(my.ATAC.seq.data.2[,c(4,5)],1,mean) + 0.01)) + 0.01 )
                                  )
rownames(my.ATAC.seq.data.3) <- my.ATAC.seq.data.3$Gene.Name

my.ATAC.seq.features <- my.ATAC.seq.data.3[intersect(my.ATAC.seq.data.3$Gene.Name, rownames(my.expressed.genes)), ]

save(my.ATAC.seq.features, 
     file = paste(Sys.Date(),"ATACseq_features_for_ML_AGING.RData", sep = "_"))


#######################
sink(file = paste(Sys.Date(),"Machine_Learning_Feature_engineering_R_session_Info_AGING.txt", sep =""))
sessionInfo()
sink()

