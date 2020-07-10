setwd('/Volumes/BB_Home//Neutrophils/Neutrophil_analysis/Lipidomics/Limma_New/')
options(stringsAsFactors = F)
library(beeswarm)

source('Functions_lipidomics_limma.R')

# Analyze Kevin's Lipidomics data
# Normalize by Shay's BCA readings
# Use protein level as covariate

# add a global data-driven normalization, which is supposed to help with robustness in metabolomics
# VSN : 
#    - https://www.nature.com/articles/srep38881
#    - https://academic.oup.com/bioinformatics/article/30/15/2155/2390365


##### NOTE: 
#  -- positive logFC means higher in males
#  -- negative logFC means higher in females

my.outprefix <- paste(Sys.Date(),"Neutrophil_Lipidomics_VSN_LION",sep="_")
my.fdr <- 0.05

#########################################################
##################  A. Preprocess data ##################  
#########################################################

####################################
# 1. read raw data and set column order
my.lipid.FA.data.1 <- read.csv('Lipidomics_Mouse_Neutrophils_BB_clean.txt', header = T, sep = "\t")
my.f.y    <- grep("F_Young", colnames(my.lipid.FA.data.1))
my.m.y    <- grep("M_Young", colnames(my.lipid.FA.data.1))
my.f.o    <- grep("F_Old",colnames(my.lipid.FA.data.1))
my.m.o    <- grep("M_Old",colnames(my.lipid.FA.data.1))

# subset and reorder
my.lipid.FA.data <- my.lipid.FA.data.1[,c(1:6,my.f.y,my.f.o,my.m.y,my.m.o)]
rownames(my.lipid.FA.data) <- my.lipid.FA.data$Lipid_ID

####################################
# 2. read BCA data and perform naive normalization to BCA
my.bca.data <- read.csv('BCA_neutrophils_Berenice_v2.txt', header = T, sep = "\t")
rownames(my.bca.data) <- my.bca.data$SampleID
my.bca.data.sorted <- my.bca.data[colnames(my.lipid.FA.data)[-c(1:6)],]

my.lipid.data.norm.bca           <- 1e6 * normalize_to_BCA(my.lipid.FA.data[-c(1:6)],my.bca.data.sorted$Total.Protein..ug.)
rownames(my.lipid.data.norm.bca) <- rownames(my.lipid.FA.data[-c(1:6)])
colnames(my.lipid.data.norm.bca) <- colnames(my.lipid.FA.data[-c(1:6)])

####################################
# 3. perform quantile/median normalization
my.lipid.data.vsn <- as.data.frame(normalizeVSN(as.matrix(my.lipid.data.norm.bca)))
rownames(my.lipid.data.vsn) <- rownames(my.lipid.FA.data[-c(1:6)])
colnames(my.lipid.data.vsn) <- colnames(my.lipid.FA.data[-c(1:6)])

# color-code 
my.colors <- rep("deeppink",16)
my.colors[grep("M_Young",colnames(my.lipid.data.vsn))] <- "deepskyblue"
my.colors[grep("M_Old",colnames(my.lipid.data.vsn))] <- "deepskyblue4"
my.colors[grep("F_Old",colnames(my.lipid.data.vsn))] <- "deeppink4"

pdf(paste(Sys.Date(),"Neutrophils_Aging_Sex_Lipidomics_Pre_PostBCA_VSN_norm_boxplots.pdf",sep = "_"), width = 15, height = 5)
par(mfrow = c(1,3))
boxplot(my.lipid.FA.data[-c(1:6)], log = 'y', las = 2, col = my.colors)
boxplot(my.lipid.data.norm.bca, log = 'y', las = 2, col = my.colors)
boxplot(my.lipid.data.vsn, las = 2, col = my.colors)
par(mfrow = c(1,1))
dev.off()

#########################################################
##################  B. Run analysis    ##################  
#########################################################

#############################################################
##%%%% 1. do MDS and PCA analysis

# MDS analysis
mds.result <- cmdscale(1-cor(my.lipid.data.vsn,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
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
my.pos.var <- apply(my.lipid.data.vsn,1,var) > 0
my.pca <- prcomp(t(my.lipid.data.vsn[my.pos.var,]),scale = TRUE)
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
dataDesign = data.frame( row.names = colnames( my.lipid.data.vsn ), 
                         age = my.Age,
                         sex = my.Sex,
                         protein = my.Protein)

# Set null and alternative models
my.model = model.matrix(~ age + sex + protein, data = dataDesign)
my.fit.all <- lmFit(my.lipid.data.vsn, design=my.model)
my.fit.all <- eBayes(my.fit.all)

my.sig.sex <- topTable(my.fit.all,coef='sexM', p.value = 1 , number = Inf)
my.sig.age <- topTable(my.fit.all,coef='age', p.value = 1 , number = Inf)

#############################################################
##%%%% a. Sex effect, age as covariate
sig.FvsM.names <- rownames(my.sig.sex[my.sig.sex$adj.P.Val < my.fdr,])

# plot significant heatmap
pdf(paste(my.outprefix,"_Heatmap_FvsM_FDR", 100* my.fdr, ".pdf", sep =""), onefile = F)
my.heatmap.title <- paste("F vs M significant (FDR<", 100* my.fdr, "%), ", length(sig.FvsM.names), " features",sep="")
pheatmap(my.lipid.data.vsn[sig.FvsM.names,], 
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

#output LION ID
my.sig.sex <- cbind(rownames(my.sig.sex),my.sig.sex)
colnames(my.sig.sex)[1] <- "Lipid_ID"

my.sig.sex.annot <- merge(my.lipid.FA.data[,1:2],my.sig.sex, by = "Lipid_ID")
my.sig.sex.annot <- my.sig.sex.annot[,c("Lipid_ID","Lipid_ID_LION","AveExpr","logFC","t","P.Value","adj.P.Val")]

write.table(my.sig.sex.annot[my.sig.sex.annot$adj.P.Val < my.fdr,], file = paste(my.outprefix,"_Limma_FvsM_LION_ID_FDR", 100* my.fdr, ".txt", sep =""), quote = F, sep = "\t", row.names = F)
write.table(my.sig.sex.annot, file = paste(my.outprefix,"_Limma_FvsM_LION_ID_ALL", ".txt", sep =""), quote = F, sep = "\t", row.names = F)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#############################################################
##%%%% b. Age effect, sex as covariate
sig.age.names <- rownames(my.sig.age[my.sig.age$adj.P.Val < my.fdr,])
sig.age.names
#### character(0)

#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()



