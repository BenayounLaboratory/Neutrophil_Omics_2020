setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/ATAC-seq_analysis/NucleoATAC/')
options(stringsAsFactors = F)
library(NucleoATACR)
library(GenomicRanges)

# Analyze NucleoATAC output
# devtools::install_github("GreenleafLab/NucleoATACR")

#######################################################
# Read in nucleosome positions:
nucs.4mF  <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_4m_F_NucleoATAC.nucmap_combined.bed.gz")
nucs.4mM  <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_4m_M_NucleoATAC.nucmap_combined.bed.gz")
nucs.21mF <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_21m_F_NucleoATAC.nucmap_combined.bed.gz")
nucs.21mM <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_21m_M_NucleoATAC.nucmap_combined.bed.gz")

get_nuc_distances <- function(nuc.calls, max_dist = 500) { # modify get_dist_between_calls to get non tabulated outuput
  sep          <- BiocGenerics::start(nuc.calls)[2:(length(nuc.calls))] - BiocGenerics::end(nuc.calls)[1:(length(nuc.calls) - 1)]
  chrom_same   <- which(as.character(GenomicRanges::seqnames(nuc.calls))[1:(length(nuc.calls) - 1)] == as.character(GenomicRanges::seqnames(nuc.calls))[2:(length(nuc.calls))])
  close.nucs   <- intersect(which(sep < max_dist), chrom_same)
  return(sep[close.nucs])
}

my.dist.nucs.4mF  <- get_nuc_distances(nucs.4mF)
my.dist.nucs.4mM  <- get_nuc_distances(nucs.4mM)
my.dist.nucs.21mF <- get_nuc_distances(nucs.21mF)
my.dist.nucs.21mM <- get_nuc_distances(nucs.21mM)

my.distances <- list("YF" = my.dist.nucs.4mF,
                     "YM" = my.dist.nucs.4mM,
                     "OF" = my.dist.nucs.21mF,
                     "OM" = my.dist.nucs.21mM)

lapply(my.distances,median)
# $YF 278
# $YM 270
# $OF 270
# $OF 265

my.test.sex.young <- wilcox.test(my.distances$YF,my.distances$YM) 
my.test.sex.old   <- wilcox.test(my.distances$OF,my.distances$OM) 
my.test.sex.young$p.value # 7.309365e-168
my.test.sex.old$p.value   #  5.898146e-14

my.test.age.f   <- wilcox.test(my.distances$YF,my.distances$OF)  
my.test.age.m   <- wilcox.test(my.distances$YM,my.distances$OM)
my.test.age.f$p.value   # 2.332269e-207
my.test.age.m$p.value   # 9.416044e-33


pdf(paste0(Sys.Date(),"_Boxplot_Dyad_distance_All_Nucs.pdf"), height = 5, width = 3.5)
boxplot(my.distances, 
        col = c("deeppink", "deepskyblue", "deeppink3","deepskyblue3"), 
        las = 1, outline = F,
        ylim = c(100,600),
        ylab = "Inter-dyad distance (bp)")
text(1.5,525, signif(my.test.sex.young$p.value,3), cex = 0.75)
text(3.5,525, signif(my.test.sex.old$p.value,3)  , cex = 0.75)
text(2,  600, signif(my.test.age.f$p.value,3)    , cex = 0.75, col = "deeppink")
text(3,  600, signif(my.test.age.m$p.value,3)    , cex = 0.75, col = "deepskyblue")
dev.off()


### Occupancy (all called nucleosomes)
my.test.sex.young <- wilcox.test(nucs.4mF$occ ,nucs.4mM$occ) 
my.test.sex.old   <- wilcox.test(nucs.21mF$occ,nucs.21mM$occ) 
my.test.sex.young$p.value # 0
my.test.sex.old$p.value   # 1.027354e-56

my.test.age.f   <- wilcox.test(nucs.4mF$occ,nucs.21mF$occ)  
my.test.age.m   <- wilcox.test(nucs.4mM$occ,nucs.21mM$occ)
my.test.age.f$p.value   #  0
my.test.age.m$p.value   # 2.972908e-131

pdf(paste0(Sys.Date(),"_Boxplot_Nucleosome_Occupancy_Scores_All_Nucs.pdf"), height = 5, width = 3.5)
boxplot(list("YF" = nucs.4mF$occ     ,
             "YM" = nucs.4mM$occ     ,
             "OF" = nucs.21mF$occ    ,
             "OM" = nucs.21mM$occ    ),
        col = c("deeppink", "deepskyblue", "deeppink3","deepskyblue3"),
        las = 1, outline = F,
        ylab = "Nucleosome Occupancy Score",
        ylim = c(0,1.2) )
text(1.5,1.1, signif(my.test.sex.young$p.value,3), cex = 0.75)
text(3.5,1.1, signif(my.test.sex.old$p.value,3)  , cex = 0.75)
text(2,  1.2, signif(my.test.age.f$p.value,3)    , cex = 0.75, col = "deeppink")
text(3,  1.2, signif(my.test.age.m$p.value,3)    , cex = 0.75, col = "deepskyblue")
dev.off() 

### the non-occ called nucleosome occupancy may be wrongly called? 
# consistent with below, 

########################################################
# read the nucleosome output with fuzziness
nucs.pos.4mF    <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_4m_F_NucleoATAC.nucpos.bed.gz")
nucs.pos.4mM    <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_4m_M_NucleoATAC.nucpos.bed.gz")
nucs.pos.21mF   <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_21m_F_NucleoATAC.nucpos.bed.gz")
nucs.pos.21mM   <- readNucs("NucleoATAC_output/Neutrophil_ATACseq_21m_M_NucleoATAC.nucpos.bed.gz")

### Fuzziness
my.test.sex.young <- wilcox.test(nucs.pos.4mF$fuzz,nucs.pos.4mM$fuzz) 
my.test.sex.old   <- wilcox.test(nucs.pos.21mF$fuzz,nucs.pos.21mM$fuzz) 
my.test.sex.young$p.value # 0
my.test.sex.old$p.value   # 7.675069e-31

my.test.age.f   <- wilcox.test(nucs.pos.4mF$fuzz,nucs.pos.21mF$fuzz)  
my.test.age.m   <- wilcox.test(nucs.pos.4mM$fuzz,nucs.pos.21mM$fuzz)
my.test.age.f$p.value   # 2.40857e-109
my.test.age.m$p.value   # 1.56614e-09

pdf(paste0(Sys.Date(),"_Boxplot_Nucleosome_Fuzziness.pdf"), height = 5, width = 3.5)
boxplot(list("YF" = nucs.pos.4mF$fuzz     ,
             "YM" = nucs.pos.4mM$fuzz     ,
             "OF" = nucs.pos.21mF$fuzz    ,
             "OM" = nucs.pos.21mM$fuzz    ),
        col = c("deeppink", "deepskyblue", "deeppink3","deepskyblue3"),
        las = 1, outline = F,
        ylab = "Nucleosome Fuzziness Score",
        ylim = c(15,55) )
text(1.5,50, signif(my.test.sex.young$p.value,3), cex = 0.75)
text(3.5,50, signif(my.test.sex.old$p.value,3)  , cex = 0.75)
text(2,  55, signif(my.test.age.f$p.value,3)    , cex = 0.75, col = "deeppink")
text(3,  55, signif(my.test.age.m$p.value,3)    , cex = 0.75, col = "deepskyblue")
dev.off() 

### Occupancy (high confidence calls only)
my.test.sex.young <- wilcox.test(nucs.pos.4mF$occ,nucs.pos.4mM$occ) 
my.test.sex.old   <- wilcox.test(nucs.pos.21mF$occ,nucs.pos.21mM$occ) 
my.test.sex.young$p.value # 7.481553e-08
my.test.sex.old$p.value   # 2.284273e-98

my.test.age.f   <- wilcox.test(nucs.pos.4mF$occ,nucs.pos.21mF$occ)  
my.test.age.m   <- wilcox.test(nucs.pos.4mM$occ,nucs.pos.21mM$occ)
my.test.age.f$p.value   #  2.892213e-58
my.test.age.m$p.value   # 5.570951e-231

pdf(paste0(Sys.Date(),"_Boxplot_Nucleosome_Occupancy_Scores_occNucs_only.pdf"), height = 5, width = 3.5)
boxplot(list("YF" = nucs.pos.4mF$occ     ,
             "YM" = nucs.pos.4mM$occ     ,
             "OF" = nucs.pos.21mF$occ    ,
             "OM" = nucs.pos.21mM$occ    ),
        col = c("deeppink", "deepskyblue", "deeppink3","deepskyblue3"),
        las = 1, outline = F,
        ylab = "Nucleosome Occupancy Score",
        ylim = c(0,1) )
text(1.5,0.9, signif(my.test.sex.young$p.value,3), cex = 0.75)
text(3.5,0.9, signif(my.test.sex.old$p.value,3)  , cex = 0.75)
text(2,  1, signif(my.test.age.f$p.value,3)    , cex = 0.75, col = "deeppink")
text(3,  1, signif(my.test.age.m$p.value,3)    , cex = 0.75, col = "deepskyblue")
dev.off() 



#######################################################
# Read in vplot and plot:
v.4mF <- read_vplot("Neutrophil_ATACseq_4m_F_NucleoATAC.VMat")
v.4mM <- read_vplot("Neutrophil_ATACseq_4m_M_NucleoATAC.VMat")
v.21mF <- read_vplot("Neutrophil_ATACseq_21m_F_NucleoATAC.VMat")
v.21mM <- read_vplot("Neutrophil_ATACseq_21m_M_NucleoATAC.VMat")

pdf(paste0(Sys.Date(),"_4mF_vplot.pdf"))
plotV(v.4mF,palette = "YlOrRd")
dev.off()

pdf(paste0(Sys.Date(),"_4mM_vplot.pdf"))
plotV(v.4mM,palette = "YlOrRd")
dev.off()

pdf(paste0(Sys.Date(),"_21mF_vplot.pdf"))
plotV(v.21mF,palette = "YlOrRd")
dev.off()

pdf(paste0(Sys.Date(),"_21mM_vplot.pdf"))
plotV(v.21mM,palette = "YlOrRd")
dev.off()
