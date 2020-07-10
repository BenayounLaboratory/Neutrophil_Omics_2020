setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Metabolomics/Limma_NEW/')
options(stringsAsFactors = F)

source('Functions_metabolomics_limma_v2.R')

# Analyze Kevin's Metabolomics data
# Normalize by Shay's BCA readings
# Use protein level as covariate

# add a global data-driven normalization, which is supposed to help with robustness:
# Quantile Normalozation: PMID: 23808607, PMID: 22878636 
# VSN : https://www.nature.com/articles/srep38881, https://academic.oup.com/bioinformatics/article/30/15/2155/2390365

##### NOTE: 
#  -- positive logFC means higher in males
#  -- negative logFC means higher in females

my.outprefix <- paste(Sys.Date(),"Neutrophil_Metabolomics_VSN",sep="_")
my.fdr <- 0.05

#########################################################
##################  A. Preprocess data ##################  
#########################################################

####################################
# 1. read raw data and set column order
my.metab.data.1 <- read.csv('BB_Mouse_Neutrophils_metabolomics.txt', header = T, sep = "\t")
my.f.y    <- grep("F_4m", colnames(my.metab.data.1))
my.m.y    <- grep("M_4m", colnames(my.metab.data.1))
my.f.o    <- grep("F_20m",colnames(my.metab.data.1))
my.m.o    <- grep("M_20m",colnames(my.metab.data.1))
colnames(my.metab.data.1)[1] <- "Compound"

# subset and reorder
my.metab.data <- my.metab.data.1[,c(1,my.f.y,my.f.o,my.m.y,my.m.o)]
rownames(my.metab.data) <- my.metab.data$Compound


####################################
# 2. read BCA data and perform naive normalization to BCA
my.bca.data <- read.csv('BCA_neutrophils_Berenice_v2.txt', header = T, sep = "\t")
rownames(my.bca.data) <- my.bca.data$SampleID_old
my.bca.data.sorted <- my.bca.data[colnames(my.metab.data)[-1],]

my.metab.data.norm.bca <- normalize_to_BCA(my.metab.data[-1],my.bca.data.sorted$Total.Protein..ug.)
rownames(my.metab.data.norm.bca) <- rownames(my.metab.data[-1])
colnames(my.metab.data.norm.bca) <- colnames(my.metab.data[-1])


####################################
# 3. perform quantile/median normalization
my.metab.data.vsn <- as.data.frame(normalizeVSN(as.matrix(my.metab.data.norm.bca)))
rownames(my.metab.data.vsn) <- rownames(my.metab.data[-1])
colnames(my.metab.data.vsn) <- colnames(my.metab.data[-1])

# color-code 
my.colors <- rep("deeppink",16)
my.colors[grep("M_4m",colnames(my.metab.data.vsn))] <- "deepskyblue"
my.colors[grep("M_20m",colnames(my.metab.data.vsn))] <- "deepskyblue4"
my.colors[grep("F_20m",colnames(my.metab.data.vsn))] <- "deeppink4"

pdf(paste(Sys.Date(),"Neutrophils_Aging_Sex_Metabolomics_Pre_PostBCA_VSN_norm_boxplots.pdf",sep = "_"), width = 15, height = 5)
par(mfrow = c(1,3))
boxplot(my.metab.data[-1], log = 'y', las = 2, col = my.colors)
boxplot(my.metab.data.norm.bca, log = 'y', las = 2, col = my.colors)
boxplot(my.metab.data.vsn, log = 'y', las = 2, col = my.colors)
par(mfrow = c(1,1))
dev.off()


####################################
# 4. extract annotatation data from matrix
my.metab.annot <- my.metab.data.1[,c("Compound","Formula","Compounds_ID", "Metabolite", "HMDB" , "KEGG")]

######
my.HMDB.split <- strsplit(my.metab.annot$HMDB,"|", fixed =T)
my.HMDBs <- vector(mode = "character", length = length(my.HMDB.split))

for (i in 1:length(my.HMDB.split)) {
  my.HMDBs[i] <- my.HMDB.split[[i]][1]
}
my.HMDBs[is.na(my.HMDBs)] <- ""
my.metab.annot$HMDB_top <- my.HMDBs

######
my.Metabolite.split <- strsplit(my.metab.annot$Metabolite,"|", fixed =T)
my.Metabolite <- vector(mode = "character", length = length(my.Metabolite.split))
for (i in 1:length(my.Metabolite.split)) {
  my.Metabolite[i] <- my.Metabolite.split[[i]][1]
}
my.Metabolite[is.na(my.Metabolite)] <- ""
my.metab.annot$Metab_top <- my.Metabolite

my.metab.annot.2 <- my.metab.annot[,c("Compound","Formula","Compounds_ID", "Metab_top", "HMDB_top")]


# retrieve annotations
my.metab.data.vsn.2 <- cbind(rownames(my.metab.data.vsn),my.metab.data.vsn)
colnames(my.metab.data.vsn.2)[1] <- "Compound"
my.metab.data.vsn.annot <- merge(my.metab.annot.2,my.metab.data.vsn.2, by = "Compound")
write.table(my.metab.data.vsn.annot, file = paste(my.outprefix,"_BCA_VSN_normalized.txt", sep =""), quote = F, sep = "\t", row.names = F)

#########################################################
##################  B. Run analysis    ##################  
#########################################################

#############################################################
##%%%% 1. do MDS and PCA analysis

# MDS analysis
mds.result <- cmdscale(1-cor(my.metab.data.vsn,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix,"MDS_plot.pdf", sep ="_")
pdf(my.mds.out)
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",
     cex=3,
     pch = 16, 
     col = my.colors,
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()

# PCA analysis
my.pos.var <- apply(my.metab.data.vsn,1,var) > 0
my.pca <- prcomp(t(my.metab.data.vsn[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix,"PCA_plot.pdf", sep ="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, 
     col= my.colors, 
     pch=16,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5) 
dev.off()


#############################################################
##%%%% 2. do limma analysis 
my.Sex  <- c(rep("F",8),rep("M",8))
my.Age  <- c(rep(4,4),rep(20,4),rep(4,4),rep(20,4))
my.Protein  <- my.bca.data.sorted$Total.Protein..ug.

# build design matrix
dataDesign = data.frame( row.names = colnames( my.metab.data.vsn ), 
                         age = my.Age,
                         sex = my.Sex,
                         protein = my.Protein)

# Set null and alternative models
my.model = model.matrix(~ age + sex + protein, data = dataDesign)
my.fit.all <- lmFit(my.metab.data.vsn, design=my.model)
my.fit.all <- eBayes(my.fit.all)

my.sig.sex <- topTable(my.fit.all,coef='sexM', p.value = 1 , number = Inf)
my.sig.age <- topTable(my.fit.all,coef='age', p.value = 1 , number = Inf)

#############################################################
##%%%% a. Sex effect, age as covariate
sig.FvsM.names <- rownames(my.sig.sex[my.sig.sex$adj.P.Val < my.fdr,])

# plot significant heatmap
pdf(paste(my.outprefix,"_Heatmap_FvsM_FDR", 100* my.fdr, ".pdf", sep =""), onefile = F)
my.heatmap.title <- paste("F vs M significant (FDR<", 100* my.fdr, "%), ", length(sig.FvsM.names), " features",sep="")
pheatmap(my.metab.data.vsn[sig.FvsM.names,], 
         scale = 'row',
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         border_color = NA,
         show_rownames = T,
         main = my.heatmap.title,
         cellwidth = 20
)
dev.off()

# output useful info to file
write.table(my.sig.sex[my.sig.sex$adj.P.Val < my.fdr,], file = paste(my.outprefix,"_Limma_FvsM_FDR", 100* my.fdr, ".txt", sep =""), quote = F, sep = "\t", row.names = T)
write.table(my.sig.sex, file = paste(my.outprefix,"_Limma_FvsM_ALL", ".txt", sep =""), quote = F, sep = "\t", row.names = T)

my.sig.sex <- cbind(rownames(my.sig.sex),my.sig.sex)
colnames(my.sig.sex)[1] <- "Compound"

my.sig.sex.annot <- merge(my.sig.sex,my.metab.annot.2, by = "Compound")

write.table(my.sig.sex.annot[my.sig.sex.annot$adj.P.Val < my.fdr,], file = paste(my.outprefix,"_Limma_FvsM_ANNOT_FDR", 100* my.fdr, ".txt", sep =""), quote = F, sep = "\t", row.names = F)
write.table(my.sig.sex.annot, file = paste(my.outprefix,"_Limma_FvsM_ANNOT_ALL", ".txt", sep =""), quote = F, sep = "\t", row.names = F)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#############################################################
##%%%% b. Age effect, sex as covariate
sig.age.names <- rownames(my.sig.age[my.sig.age$adj.P.Val < my.fdr,])
#"8.04_244.1112m/z"
#                      logFC  AveExpr        t      P.Value  adj.P.Val        B
# 8.04_244.1112m/z 0.1101104 4.905809 6.353475 1.567115e-05 0.03847266 3.305435

my.sig.age <- cbind(rownames(my.sig.age),my.sig.age)
colnames(my.sig.age)[1] <- "Compound"

my.sig.age.annot <- merge(my.sig.age,my.metab.annot.2, by = "Compound")

my.sig.age.annot[my.sig.age.annot$adj.P.Val < my.fdr,]
#          Compound     logFCAveExpr   t         P.Value   adj.P.Val        B     Formula  Compounds_ID    Metab_top  HMDB_top
#2135 8.04_244.1112m/z 0.1101104 4.905809 6.353475 1.567115e-05 0.03847266 3.305435 C10H17N3O2S  2135 pHILIC_244.1112_8 Biotin_amide HMDB01458
#############################################################



#############################################################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()

