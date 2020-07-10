setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Age_Expression/Features/Feature_Engineering')
options(stringsAsFactors = F)

## run similar to the sex ML
load('2020-06-01_RNA_features_for_ML_AGING.RData')
load('2020-06-01_ATACseq_features_for_ML_AGING.RData')
load('2020-06-01_summary_TF_features_from_GMT_for_ML_AGING.RData')
load('2020-06-01_HOMER_sequence_features_for_ML_AGING.RData')


# Add gene name column for RNA and TF features
my.RNA.features$Gene.Name         <- rownames(my.RNA.features)
my.TF.features.FINAL$Gene.Name    <- rownames(my.TF.features.FINAL)

rownames(my.ATAC.seq.features)    <- my.ATAC.seq.features$Gene.Name
rownames(my.seq.feat.alls)        <- my.seq.feat.alls$Gene.Name

# Merge dataframe sequentially
my.RNA.ATAC.features             <- merge(my.RNA.features              , my.ATAC.seq.features    , by = "Gene.Name")
my.RNA.ATAC.seq.features         <- merge(my.RNA.ATAC.features         , my.seq.feat.alls        , by = "Gene.Name")
my.RNA.ATAC.seq.TF.features      <- merge(my.RNA.ATAC.seq.features     , my.TF.features.FINAL    , by = "Gene.Name")

summary(my.RNA.ATAC.features$AGING)
# DOWN   UP 
# 1840 1581 

#############################################
### rename and save
# remove "FDR","logFC","Chrosomome" features for learning (keep chromosome info as autosome vs Sex chromome only)
my.ALL.features.2group                <- my.RNA.ATAC.seq.TF.features[,!(colnames(my.RNA.ATAC.seq.TF.features) %in% c("FDR","logFC","Chrosomome"))]

summary(my.ALL.features.2group$AGING)
# DOWN   UP 
# 1840 1581 

save(my.ALL.features.2group,
     file = paste0(Sys.Date(),"_ALL_ML_features_for_AGING_Expression.RData"))
         
         
