setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/ATAC-seq_analysis/NucleoATAC/')
options(stringsAsFactors = F)

# 2020-05-01
# Analyze fragment distribution at ATAC_peaks

#######################################################
# Read in nucleosome positions:
Peak.fragsizes  <- read.table("Neutrophil_ATACseq_Peaks_meta.QCnoMt.FragSizes.ANNOTATED.txt", header = T, sep = "\t")
colnames(Peak.fragsizes)

Peak.fragsizes.median <- Peak.fragsizes[,c(1:7,grep("_MedFragSize",colnames(Peak.fragsizes)))]

# get median value over replicates
my.distances <- list("YF" = apply(Peak.fragsizes.median[,grep("4m_F", colnames(Peak.fragsizes.median))],1,median) ,
                     "YM" = apply(Peak.fragsizes.median[,grep("4m_M", colnames(Peak.fragsizes.median))],1,median) ,
                     "OF" = apply(Peak.fragsizes.median[,grep("21m_F",colnames(Peak.fragsizes.median))],1,median) ,
                     "OM" = apply(Peak.fragsizes.median[,grep("21m_M",colnames(Peak.fragsizes.median))],1,median) )

median(my.distances$YF, na.rm = T) #  193.5
median(my.distances$YM, na.rm = T) #  182.5
median(my.distances$OF, na.rm = T) #  190
median(my.distances$OM, na.rm = T) #  181.5

my.test.sex.young <- wilcox.test(my.distances$YF, my.distances$YM) 
my.test.sex.old   <- wilcox.test(my.distances$OF,my.distances$OM) 
my.test.sex.young$p.value #  1.493039e-51
my.test.sex.old$p.value   #  1.190874e-48

my.test.age.f   <- wilcox.test(my.distances$YF,my.distances$OF)  
my.test.age.m   <- wilcox.test(my.distances$YM,my.distances$OM)
my.test.age.f$p.value # 1.694868e-44
my.test.age.m$p.value # 1.03191e-31

pdf(paste0(Sys.Date(),"_Peaks_Median_FragSizes.pdf"), height = 5, width = 3.5)
boxplot(my.distances, 
        col = c("deeppink", "deepskyblue", "deeppink3","deepskyblue3"), 
        las = 1, outline = F,
        ylim = c(0,700),
        ylab = "Median Fragment size (bp)")
text(1.5,600, signif(my.test.sex.young$p.value,3), cex = 0.75)
text(3.5,600, signif(my.test.sex.old$p.value,3)  , cex = 0.75)
text(2,  700, signif(my.test.age.f$p.value,3)    , cex = 0.75, col = "deeppink")
text(3,  700, signif(my.test.age.m$p.value,3)    , cex = 0.75, col = "deepskyblue")
dev.off()



