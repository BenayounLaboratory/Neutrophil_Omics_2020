setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/GSEA/Bubble_Chart/')
options(stringsAsFactors = F)

source('Plot_bubble_chart_functions_v4.R')

# try to get a divergent bubble plot going (top each side)
# using the GSEA Data
# max top 10 most significant for each sex

########################################################################################################################
### SEX
get_sex_bubble_plot("Neutrophils_Sex_MSigDB_KEGG"                      , '../GSEA_Sex/2020-05-20_Neutrophil_Sex_Effect_MSigDB_KEGG_FDR_5_Phenotest_GSEA_Analysis_table_79_significant.txt'            )
get_sex_bubble_plot("Neutrophils_Sex_MSigDB_Reactome"                  , '../GSEA_Sex/2020-05-20_Neutrophil_Sex_Effect_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_643_significant.txt'       )
get_sex_bubble_plot("Neutrophils_Sex_ENS.GO.ALL"                       , '../GSEA_Sex/2020-05-20_Neutrophil_Sex_Effect_ENS.GO.ALL_FDR_5_Phenotest_GSEA_Analysis_table_2118_significant.txt'           )
get_sex_bubble_plot("Neutrophils_Sex_TF_targets"                       , '../GSEA_Sex/2020-05-20_Neutrophil_Sex_Effect_TF_Targets_ALL_FDR_5_Phenotest_GSEA_Analysis_table_71_significant.txt'        )


########################################################################################################################
#### AGE
get_age_bubble_plot("Neutrophils_Aging_MSigDB_KEGG"                      , '../GSEA_Age/2020-05-20_Neutrophil_Aging_Effect_MSigDB_KEGG_FDR_5_Phenotest_GSEA_Analysis_table_28_significant.txt'        )
get_age_bubble_plot("Neutrophils_Aging_MSigDB_Reactome"                  , '../GSEA_Age/2020-05-20_Neutrophil_Aging_Effect_MSigDB_Reactome_FDR_5_Phenotest_GSEA_Analysis_table_446_significant.txt'   )
get_age_bubble_plot("Neutrophils_Aging_ENS.GO.ALL"                       , '../GSEA_Age/2020-05-20_Neutrophil_Aging_Effect_ENS.GO.ALL_FDR_5_Phenotest_GSEA_Analysis_table_1718_significant.txt'       )
get_age_bubble_plot("Neutrophils_Aging_TF_targets"                       , '../GSEA_Age/2020-05-20_Neutrophil_Aging_Effect_TF_Targets_ALL_FDR_5_Phenotest_GSEA_Analysis_table_102_significant.txt'    )
