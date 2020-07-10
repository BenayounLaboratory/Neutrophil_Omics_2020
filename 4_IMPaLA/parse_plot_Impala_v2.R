setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Integration/Impala/30-06_validated_metab_lipids/')
options(stringsAsFactors = F)

library(bitops)

# 2020-06-30
# parse intergrated analysis with Impala
# validated metabolites only
# Lipids with HMDB ID

my.impala.male   <- read.csv('2020-06-30_Male_genes_metab_lipid_ORA_results.csv', header = T)
my.impala.female <- read.csv('2020-06-30_Female_genes_metab_lipid_ORA_results.csv', header = T)

# Only retain pathways that are significant overall at FDR 5%
my.impala.male.overall   <- my.impala.male[my.impala.male$Q_joint < 0.05,]     # 512
my.impala.female.overall <- my.impala.female[my.impala.female$Q_joint < 0.05,] # 173
  
# Only retain pathways that are significant for metabolites, for genes, and jointly
my.impala.male.filt   <- my.impala.male.overall[ bitAnd(my.impala.male.overall$Q_genes < 0.05, my.impala.male.overall$Q_metabolites < 0.05)>0, ]
my.impala.female.filt <- my.impala.female.overall[ bitAnd(my.impala.female.overall$Q_genes < 0.05, my.impala.female.overall$Q_metabolites < 0.05)>0, ]

nrow(my.impala.male.filt)    # 33
nrow(my.impala.female.filt)  # 1

write.table(my.impala.male.filt, file = paste0(Sys.Date(),"_Impala_FDR5_overall_genes_metabs_MALES.txt"), sep = "\t", quote = F, row.names = F)
write.table(my.impala.female.filt, file = paste0(Sys.Date(),"_Impala_FDR5_overall_genes_metabs_FEMALES.txt"), sep = "\t", quote = F, row.names = F)

write.table(my.impala.male.overall, file = paste0(Sys.Date(),"_Impala_FDR5_overall_MALES.txt"), sep = "\t", quote = F, row.names = F)
write.table(my.impala.female.overall, file = paste0(Sys.Date(),"_Impala_FDR5_overall_FEMALES.txt"), sep = "\t", quote = F, row.names = F)

# create a variable with info on whether pathway is in both or only one "omic" realm
my.impala.male.overall$Q_BOTH   <- "Genes_Only"
my.impala.female.overall$Q_BOTH <- "Genes_Only"

my.impala.male.overall$Q_BOTH[my.impala.male.overall$Q_metabolites < 0.05]     <- "Metabolites_Only"
my.impala.female.overall$Q_BOTH[my.impala.female.overall$Q_metabolites < 0.05] <- "Metabolites_Only"

my.impala.male.overall$Q_BOTH[ bitAnd(my.impala.male.overall$Q_genes < 0.05, my.impala.male.overall$Q_metabolites < 0.05)>0 ]       <- "Both"
my.impala.female.overall$Q_BOTH[ bitAnd(my.impala.female.overall$Q_genes < 0.05, my.impala.female.overall$Q_metabolites < 0.05)>0 ] <- "Both"

my.impala.male.overall$Sex_bias    <- "Male"
my.impala.female.overall$Sex_bias  <- "Female"

my.impala.summary <- rbind(my.impala.male.overall,my.impala.female.overall)

write.table(my.impala.summary, file = paste0(Sys.Date(),"_Impala_FDR5_overall_summary.txt"), sep = "\t", quote = F, row.names = F)


