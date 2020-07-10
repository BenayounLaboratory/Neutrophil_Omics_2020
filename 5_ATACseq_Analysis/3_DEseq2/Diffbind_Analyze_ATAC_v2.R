setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/ATAC-seq_analysis/Diffbind')
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library('DiffBind')
library("rtracklayer")

# DS bam files
# New clean HOMER peaks (meta for each group)

################################################################################
####################   Diffbind Analysis Ahr female vs. Male  ####################
################################################################################
ntph.aging <- dba(sampleSheet="ATAC_Sex_Aging_Neutrophils_samples_DS.csv",skipLines=0,attributes=c(DBA_ID,DBA_CONDITION))
ntph.aging <- dba.count(ntph.aging)

pdf(paste0(Sys.Date(),"ATAC_Sex_Aging_Neutrophils_heatmap_DS.pdf"))
plot(ntph.aging, colScheme="Reds")
dev.off()

ntph.aging <- dba.contrast(ntph.aging, categories=DBA_CONDITION,minMembers=2)
ntph.aging <- dba.analyze(ntph.aging,method=DBA_ALL_METHODS)
ntph.aging

# 20 Samples, 95926 sites in matrix:
#   ID Condition Replicate Caller Intervals FRiP
# 1   4m_F_101        YF         1 counts     95926 0.53
# 2   4m_F_102        YF         2 counts     95926 0.30
# 3   4m_F_103        YF         3 counts     95926 0.74
# 4   4m_F_104        YF         4 counts     95926 0.68
# 5   4m_F_105        YF         5 counts     95926 0.59
# 6   4m_M_111        YM         1 counts     95926 0.31
# 7   4m_M_112        YM         2 counts     95926 0.54
# 8   4m_M_113        YM         3 counts     95926 0.67
# 9   4m_M_114        YM         4 counts     95926 0.69
# 10  4m_M_115        YM         5 counts     95926 0.54
# 11 21m_F_106        OF         1 counts     95926 0.80
# 12 21m_F_107        OF         2 counts     95926 0.62
# 13 21m_F_108        OF         3 counts     95926 0.32
# 14 21m_F_109        OF         4 counts     95926 0.58
# 15 21m_F_110        OF         5 counts     95926 0.22
# 16 21m_M_116        OM         1 counts     95926 0.70
# 17 21m_M_117        OM         2 counts     95926 0.27
# 18 21m_M_118        OM         3 counts     95926 0.63
# 19 21m_M_119        OM         4 counts     95926 0.55
# 20 21m_M_120        OM         5 counts     95926 0.44



my.norm.counts <- dba.peakset(ntph.aging,bRetrieve=TRUE)

write.table(data.frame(my.norm.counts), 
            file = paste0(Sys.Date(),"_ATAC_Sex_Aging_Neutrophils_Normalized_count_matrix.txt"),
            quote = F, col.names = T, row.names = F, sep = "\t")

write.table(data.frame(my.norm.counts)[,c(1:3)], 
            file = paste0(Sys.Date(),"_ATAC_Sex_Aging_Neutrophils_diffbind_peaks.bed"),
            quote = F, col.names = F, row.names = F, sep = "\t")
