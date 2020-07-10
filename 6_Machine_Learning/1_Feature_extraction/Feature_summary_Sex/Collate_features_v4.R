setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/Features/Feature_engineering')
options(stringsAsFactors = F)

load('2020-05-20_RNA_features_for_ML.RData')
load('2020-05-20_ATACseq_features_for_ML.RData')
load('2020-05-20_summary_TF_features_from_GMT_for_ML.RData')
load('2020-05-20_HOMER_sequence_features_for_ML.RData')


# Add gene name column for RNA and TF features
my.RNA.features$Gene.Name         <- rownames(my.RNA.features)
my.TF.features.FINAL$Gene.Name    <- rownames(my.TF.features.FINAL)

rownames(my.ATAC.seq.features)    <- my.ATAC.seq.features$Gene.Name
rownames(my.seq.feat.alls)        <- my.seq.feat.alls$Gene.Name
# 
# # my.ATAC.seq.features$Gene.Name
# # my.seq.feat.alls$Gene.Name

# Merge dataframe sequentially
my.RNA.ATAC.features             <- merge(my.RNA.features              , my.ATAC.seq.features    , by = "Gene.Name")
my.RNA.ATAC.seq.features         <- merge(my.RNA.ATAC.features         , my.seq.feat.alls        , by = "Gene.Name")
my.RNA.ATAC.seq.TF.features      <- merge(my.RNA.ATAC.seq.features     , my.TF.features.FINAL    , by = "Gene.Name")

summary(my.RNA.ATAC.features$SEX_DIMORPHISM)
# FEMALE   MALE 
# 1013    623 
# summary(my.RNA.ATAC.seq.features$SEX_DIMORPHISM)
# sequence features are missing for some

#############################################
### rename and save
# remove "FDR","logFC","Chrosomome" features for learning (keep chromosome info as autosome vs Sex chromome only)
my.ALL.features.2group                <- my.RNA.ATAC.seq.TF.features[,!(colnames(my.RNA.ATAC.seq.TF.features) %in% c("FDR","logFC","Chrosomome"))]

summary(my.ALL.features.2group$SEX_DIMORPHISM)
# FEMALE   MALE 
# 1013      623 

save(my.ALL.features.2group,
     file = paste0(Sys.Date(),"_ALL_ML_features_for_SexDim_Expression.RData"))
         
         
