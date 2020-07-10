setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/ATAC-seq_analysis/NucleoATAC/')
options(stringsAsFactors = F)

library('pheatmap')

# Analyze NucleoATAC output
# use HOMER bedgraph signal extraction


##########################################################################################
# read sorted data matrix
my.occ.mat <- read.table('HOMER_Neutrophil_NucleoATAC_occupancy_4mFM_21mFM.txt',sep="\t",header=T)

my.4mF.cols  <-   2:102
my.4mM.cols  <- 103:203
my.21mF.cols <- 204:304
my.21mM.cols <- 305:405


my.occ.4mF  <- my.occ.mat[,my.4mF.cols ]
my.occ.4mM  <- my.occ.mat[,my.4mM.cols ]
my.occ.21mF <- my.occ.mat[,my.21mF.cols]
my.occ.21mM <- my.occ.mat[,my.21mM.cols]

my.av.occ.4mF  <- apply(my.occ.4mF  , 2, mean)
my.av.occ.4mM  <- apply(my.occ.4mM  , 2, mean)
my.av.occ.21mF <- apply(my.occ.21mF , 2, mean)
my.av.occ.21mM <- apply(my.occ.21mM , 2, mean)

my.median.occ.4mF  <- apply(my.occ.4mF  , 2, median)
my.median.occ.4mM  <- apply(my.occ.4mM  , 2, median)
my.median.occ.21mF <- apply(my.occ.21mF , 2, median)
my.median.occ.21mM <- apply(my.occ.21mM , 2, median)

my.dist <- seq(-500,500,10)


#### Sex effect on nucleosome occupancy
ks.test(my.median.occ.4mF,  my.median.occ.4mM) # p-value = 0.0004839
ks.test(my.median.occ.21mF, my.median.occ.21mM) # p-value = 0.006672

pdf(paste0(Sys.Date(),"_Median_NucleoATAC_occupancy_at_expressed_genes_Sex.pdf"), height = 7, width = 10)
par(mfrow=c(1,2))
plot(my.dist,my.median.occ.4mF, type = 'l', 
     xlim = c(-500,500), ylim = c(0,0.55), 
     col  = "deeppink", las = 1,
     xlab = "Distance to TSS (bp)",
     ylab = "Nucleosome Occupancy Score" )
lines(my.dist,my.median.occ.4mM, col  = "deepskyblue")
text(-400,0.5,"p=0.0004839", cex = 0.75)
legend("bottomright",col=c("deeppink","deepskyblue"),c("YF","YM"),pch = "_",bty='n')
plot(my.dist,my.median.occ.21mF, type = 'l', 
     xlim = c(-500,500), ylim = c(0,0.55), 
     col  = "deeppink3", las = 1,
     xlab = "Distance to TSS (bp)",
     ylab = "Nucleosome Occupancy Score" )
lines(my.dist,my.median.occ.21mM, col  = "deepskyblue3")
text(-400,0.5,"p=0.006672", cex = 0.75)
legend("bottomright",col=c("deeppink3","deepskyblue3"),c("OF","OM"),pch = "_",bty='n')
dev.off()


ks.test(my.median.occ.4mF,  my.median.occ.21mF) # p-value = 0.0002698
ks.test(my.median.occ.4mM, my.median.occ.21mM) # p-value = 0.006672

pdf(paste0(Sys.Date(),"_Median_NucleoATAC_occupancy_at_expressed_genes_AGING.pdf"), height = 7, width = 10)
par(mfrow=c(1,2))
plot(my.dist,my.median.occ.4mF, type = 'l', 
     xlim = c(-500,500), ylim = c(0,0.55), 
     col  = "deeppink", las = 1,
     xlab = "Distance to TSS (bp)",
     ylab = "Nucleosome Occupancy Score" )
lines(my.dist,my.median.occ.21mF, col  = "deeppink3")
text(-400,0.5,"p=0.0002698", cex = 0.75)
legend("bottomright",col=c("deeppink","deeppink3"),c("YF","OF"),pch = "_",bty='n')
plot(my.dist,my.median.occ.4mM, type = 'l', 
     xlim = c(-500,500), ylim = c(0,0.55), 
     col  = "deepskyblue", las = 1,
     xlab = "Distance to TSS (bp)",
     ylab = "Nucleosome Occupancy Score" )
lines(my.dist,my.median.occ.21mM, col  = "deepskyblue3")
text(-400,0.5,"p=0.006672", cex = 0.75)
legend("bottomright",col=c("deepskyblue","deepskyblue3"),c("YM","OM"),pch = "_",bty='n')
dev.off()


