setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Metabolomics/Metaboanalyst_R')
options(stringsAsFactors = F)

library(bitops)
library(MetaboAnalystR)

# try to run metaboanalyst R using M/z results in mixed mode

#####
mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
SetPeakFormat("mpt")
mSet <- UpdateInstrumentParameters(mSet, 5, "mixed", "yes");


mSet <- Read.PeakListData(mSet, "2020-03-26_Neutrophil_metabolomics_MIXED_MODE_Males_neg.txt");
mSet <- SanityCheckMummichogData(mSet)
mSet <- SetPeakEnrichMethod(mSet, "gsea", "v2")
mSet <- SetMummichogPval(mSet, 0.05)

set.seed(123456789)
mSet <- PerformPSEA(mSet, "mmu_kegg", "current", 1000)

mSet <- PlotPeaks2Paths(mSet, paste0(Sys.Date(),"Neutrophil_sex_metaboanalyst_peaks_to_paths_GSEA_KEGG"), "pdf", 300, width=20)

my.psea.results <- data.frame(mSet$mummi.gsea.resmat)
write.table(my.psea.results, file = paste0(Sys.Date(),"Neutrophil_sex_metaboanalyst_peaks_to_paths_GSEA_KEGG_table.txt"), sep = "\t", quote= F  )
