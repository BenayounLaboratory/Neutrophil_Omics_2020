setwd("e://Metabolomics2/181010_Mouse_neutrophils_aging_BB/")

sample <- read.csv("PQI_20181010_Mouse_neutrophils_aging_RPLC_pos_new.csv")

###############################################################################################################################################
############# Remove background #############
###############################################################################################################################################
##  Create compound_ID column
sample$m.z <- round(sample$m.z, digits = 4)
sample$Retention.time..min. <- round(sample$Retention.time..min., digits = 1)
sample$Compounds_ID <- paste(sample$m.z,sample$Retention.time..min.,sep="_")

## Remove background in each batch separately. Signal must be 2 fold higher in samples compared to blks
# if missing value it gives 1 to avoid errors when the ratio is calculated
colnames(sample)
sample[,11:34] <- sample[,11:34]+1

# Batch1
sample$Mean_blk <- apply(sample[,grep("blk", colnames(sample))],1,mean)
sample$Mean_bio <- apply(sample[,13:28],1,mean)
sample$ratio <- sample$Mean_bio/sample$Mean_blk
sample <- sample[which(sample$ratio >= 2.0),]

###############################################################################################################################################
############# Clean-up of the dataset: missing values and imputation #############
###############################################################################################################################################

##  Remove unecessary samples
sample3 <- sample[,-grep("blk", colnames(sample))]

##  Remove metabolic features if missing values as follow:
# Remove metabolic features with at least one missing values in 1/6 of the QCs
#length(which(rowSums(sample3[grep("QC[[:digit:]]", colnames(sample3))] == 1) >= 1))
#sample4 <- sample3[-which(rowSums(sample3[grep("QC[[:digit:]]", colnames(sample3))] == 1) >= 1),]
sample4 <- sample3
# Remove metabolic features with missing values in more than 1/3 of the samples (here: 57 samples)
length(which(rowSums(sample4[11:26] == 1) >= 5))
sample.mval <- sample4[-which(rowSums(sample4[11:26] == 1) >= 5),]

###############################################################################################################################################
############# Correlation dilution QC #############
###############################################################################################################################################

sample_dln_filt <- sample.mval

Med_QC_dln1 <- c()
Med_QC_dln2 <- c()
Med_QC_dln4 <- c()
Med_QC_dln8 <- c()
for (i in 1:dim(sample_dln_filt)[1]) {Med_QC_dln1[i] <- sample_dln_filt[i,grep("QC1$", colnames(sample_dln_filt))]}
for (i in 1:dim(sample_dln_filt)[1]) {Med_QC_dln2[i] <- sample_dln_filt[i,grep("QC_dln2", colnames(sample_dln_filt))]}
for (i in 1:dim(sample_dln_filt)[1]) {Med_QC_dln4[i] <- sample_dln_filt[i,grep("QC_dln4", colnames(sample_dln_filt))]}
for (i in 1:dim(sample_dln_filt)[1]) {Med_QC_dln8[i] <- sample_dln_filt[i,grep("QC_dln8", colnames(sample_dln_filt))]}
cor_matrix <- data.frame(Med_QC_dln1,Med_QC_dln2,Med_QC_dln4,Med_QC_dln8)
dln <- c(1,0.5,0.25,0.125)
for (i in 1:dim(cor_matrix)[1]) {
  cor_matrix$correlation[i] <- cor(t(cor_matrix[i,1:4]),dln)^2
}
sample_dln_filt$correlation <- cor_matrix$correlation
sample_dln_filt <- sample_dln_filt[which(sample_dln_filt$correlation > 0.6),]
plot(density(cor_matrix[complete.cases(cor_matrix),]$correlation))
sample_dln_filt <- sample_dln_filt[,-grep("QC_dln", colnames(sample_dln_filt))]

rm(list = c("cor_matrix","dln","Med_QC_dln1","Med_QC_dln2","Med_QC_dln4","Med_QC_dln8"))

###############################################################################################################################################
############# Imputation #############
###############################################################################################################################################
tmp <- sample_dln_filt[,11:29]
tmp[tmp == 1] <- NA
tmp <- log2(tmp)

pdf("NA_imputation_1.8sd.pdf")
for (i in 1:ncol(tmp)){
  y <- which(is.na(tmp[,i]) == T)
  tmp[which(is.na(tmp[,i]) == T),i] <- rnorm(n = length(which(is.na(tmp[,i]) == T)), 
                                             mean = mean(as.numeric(tmp[,i]),na.rm = T)-(1.8*sd(tmp[,i],na.rm = T)), 
                                             sd = sd(tmp[,i],na.rm = T)*0.3)
  ggdata <- data.frame(Log2_Int = tmp[,i])
  ggdata$Condition <- '1'
  ggdata$Condition[y] <- '0'
  p <- ggplot(ggdata, aes(x = Log2_Int, fill = Condition)) + geom_bar(color="black",binwidth = 0.5) + 
    theme_bw() + ggtitle(colnames(tmp)[i]) + scale_fill_manual(values = c("orange", "grey"))
  print(p)
}
dev.off()

sample_dln_filt[,11:29] <- 2^tmp

###############################################################################################################################################
############# Database Identification #############
###############################################################################################################################################
sample_dln_filt$RT <- sample_dln_filt$Retention.time..min.
sample_dln_filt$mz <- sample_dln_filt$m.z
sample_dln_filt$Neutral.mass <- sample_dln_filt$mz-1.00728

library('bitops')

DB <- read.csv("e:/Metabolomics2/Database_Hannes.csv", header = TRUE)
sample_dln_filt$Metabolite <- 0
sample_dln_filt$KEGG <- 0
sample_dln_filt$HMDB <- 0
sample_dln_filt$CAS <- 0
sample_dln_filt$ChEBI <- 0
sample_dln_filt$Formula <- 0
sample_dln_filt$ppm <- 0

for (i in 1:dim(sample_dln_filt)[1]) { #for every row in sample_dln_filt file
  range <- DB[which(DB$mass < sample_dln_filt$Neutral.mass[i]+0.01 & DB$mass > sample_dln_filt$Neutral.mass[i]-0.01),] #range = subset of DB where mass matches sample_dln_filt neutral mass within +/- 0.01
  my.cur.ppms <- rep(0,dim(range)[1]) #my.cur.ppms is a n-length vector of 0s, where n is num of rows in "range".
  for (j in 1:dim(range)[1]) { # for every row in range
    my.cur.ppms[j] <- 1e6*((sample_dln_filt[i,]$Neutral.mass-range[j,]$mass)/sample_dln_filt[i,]$Neutral.mass) #my.cur.ppm for that row is 1e6 x 
  }
  my.same.mass <- which(bitAnd(my.cur.ppms < 5, my.cur.ppms > -5) > 0) # select things with ppm +/- 5 of database mass
  my.ID <- range[my.same.mass,2]
  my.KEGG <- range[my.same.mass,13]
  my.HMDB <- range[my.same.mass,7]
  my.CAS <- range[my.same.mass,12]
  my.ChEBI <- range[my.same.mass,4] 
  my.Formula <- as.character(range[my.same.mass,3])
  my.ppm <-  my.cur.ppms[my.same.mass]
  
  if(length(my.same.mass) >= 1){
    sample_dln_filt$Metabolite[i] <- paste(my.ID, collapse = "|")
    sample_dln_filt$KEGG[i] <- paste(my.KEGG, collapse = "|")
    sample_dln_filt$HMDB[i] <- paste(my.HMDB, collapse = "|")
    sample_dln_filt$CAS[i] <- paste(my.CAS, collapse = "|")
    sample_dln_filt$ChEBI[i] <- paste(my.ChEBI, collapse = "|")
    if (length(unique(my.Formula)) > 1 ){
      sample_dln_filt$Formula[i] <- paste(my.Formula, collapse = "|")
      sample_dln_filt$ppm[i] <- paste(my.ppm, collapse = "|")
    }else{
      sample_dln_filt$Formula[i] <- my.Formula[1]
      sample_dln_filt$ppm[i] <- my.ppm[1]
    }
  }else{
  }
}
length(which(sample_dln_filt$Metabolite != 0))
write.csv(sample_dln_filt,file="BB_Mouse_neutrophils_RPLC_pos.csv")

rm(list = c("DB", "my.cur.ppms", "my.same.mass", "my.ID", "my.KEGG", "my.HMDB", "my.CAS", "my.ChEBI", "my.Formula", "my.ppm"))

###############################################################################################################################################
############# Sanity checks #############
###############################################################################################################################################

##  PCA plot Batch effect
data_nz <- sample_dln_filt[,11:29] # No normalization

# Log transformation
data_nz.log10 <- log10(data_nz)

# Range scaling
data_nz.mean <- apply(data_nz.log10,1,mean)
data_nz.log10.centered <- data_nz.log10 - data_nz.mean
data_nz.max <- apply(data_nz.log10.centered,1,max)
data_nz.min <- apply(data_nz.log10.centered,1,min)
data_nz.range <- data_nz.max - data_nz.min
data_nz.norm <- data_nz.log10.centered/data_nz.range
data_nz.norm <- data_nz.norm[complete.cases(data_nz.norm),]

# Clustering using hcluster from amap package, use spearman and complete 
library('amap')
hc <- hcluster(t(data_nz.norm), method = "spearman", diag = FALSE, upper = FALSE,
               link = "complete", members = NULL, nbproc = 2,
               doubleprecision = TRUE)
plot(hc, hang = -1)

rm(list = c("df", "med1", "temp", "temp_V1","temp_V2","temp_V3","temp_V4", "data_nz", "data_nz.log10", "data_nz.mean", "data_nz.log10.centered", "data_nz.max", "data_nz.min", "data_nz.range", "data_nz.norm","data.pca","sum.pca","imp.pca","var.pca","cum.pca","cls","pc.num","pclabels","Comp","ggdata"))


