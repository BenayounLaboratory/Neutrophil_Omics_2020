setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Integration/Impala/Lists')
options(stringsAsFactors = F)

# 2020-06-30
# generate list of HMDB ids for significant lipids

my.lipid.annot <- read.csv('Kevin_Lipid_annotation.txt', sep = "\t", header = T)

my.sig.lipids <- read.csv("../../../Lipidomics/Limma_New/2020-03-23_Neutrophil_Lipidomics_VSN_LION_Limma_FvsM_LION_ID_FDR5.txt", sep = "\t", header = T)

my.sig.annot <- merge(my.lipid.annot[,c("Lipid_Name_LION","HMDB_ID")],my.sig.lipids, by.x = "Lipid_Name_LION", by.y = "Lipid_ID_LION")

my.male.lipid.hmdb   <- unique(my.sig.annot$HMDB_ID[my.sig.annot$logFC > 0])
my.female.lipid.hmdb <- unique(my.sig.annot$HMDB_ID[my.sig.annot$logFC < 0])

write.table(my.male.lipid.hmdb,   "Male_FDR5_lipids.txt", col.names = F, row.names = F, quote = F)
write.table(my.female.lipid.hmdb, "Female_FDR5_lipids.txt", col.names = F, row.names = F, quote = F)