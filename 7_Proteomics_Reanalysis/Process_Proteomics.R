setwd('/Volumes/BB_Home/Neutrophils/Public_data/Human_Neutrophil_proteomics')
options(stringsAsFactors = F)

library(readxl)
library(bitops)

# 2020-06-18
# Analyze DIA proteomics data from paper

# Read supplemental table with DIA data
my.prots <- read_xlsx("141289_1_supp_251888_pjqfps.xlsx", sheet = 2, col_names = TRUE, skip = 1)

# Keep healthy samples only
my.prots.healthy <- my.prots[my.prots$Group == "HD",]
my.prots.healthy$SampleID <- paste0(my.prots.healthy$Group,my.prots.healthy$Sample)

my.uniq.id   <- unique(my.prots.healthy$SampleID)
my.uniq.prot <- unique(my.prots.healthy$Gene)

# extract annotation columns
my.annots <- unique(my.prots.healthy[,2:4])

## parse data into DIA expression matrix
my.mat.prots <- data.frame(matrix(0,length(my.uniq.prot),length(my.uniq.id)))
rownames(my.mat.prots) <- my.uniq.prot
colnames(my.mat.prots) <- my.uniq.id

for ( i in 1:length(my.uniq.prot)) {
  for (j  in 1:length(my.uniq.id)) {
    my.prot.lines <- which(bitAnd(my.prots.healthy$Gene %in% my.uniq.prot[i],my.prots.healthy$SampleID %in% my.uniq.id[j])>0)
    my.mat.prots[i,j] <- my.prots.healthy$`Protein quantity`[my.prot.lines]
  }
}

save(my.mat.prots, my.annots, file = paste0(Sys.Date(),"_parsed_NTPH_DIA_proteomics.RData"))
##########################################################################################

##########################################################################################
# DDX3Yshould be expressed only in males
#    => use to look at correlation between sex and granule expression

load('2020-06-18_parsed_NTPH_DIA_proteomics.RData')

#################################
##### A. ELANE (ELNE in the dataset)

ddx3y.test <- cor.test(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["ELNE",]), method = "spearman") 

pdf(paste0(Sys.Date(),"_DIA_proteomics_DDX3Y_ELANE.pdf"))
plot(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["ELNE",]),
     xlab = "DIA DDX3Y levels (A.U.)",
     ylab = "DIA ELANE levels (A.U.)",
     pch = 16,
     col = "tomato",
     xlim = c(0,2e5),
     ylim = c(0,5e7),
     las = 1)
text(2e5, 4.75e7, paste0("Rho = ",signif(ddx3y.test$estimate,3)), pos = 2)
text(2e5, 4.3e7, paste0("p = ",signif(ddx3y.test$p.value,3)), pos = 2)
dev.off()

#################################
##### B. MPO (PERM in the dataset)
### PERM	Myeloperoxidase, PRTN3	Myeloblastin

ddx3y.test <- cor.test(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["PERM",]), method = "spearman") 

pdf(paste0(Sys.Date(),"_DIA_proteomics_DDX3Y_MPO.pdf"))
plot(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["PERM",]),
     xlab = "DIA DDX3Y levels (A.U.)",
     ylab = "DIA MPO levels (A.U.)",
     pch = 16,
     col = "tomato",
     xlim = c(0,2e5),
     ylim = c(0,8e7),
     las = 1)
text(2e5, 7.6e7, paste0("Rho = ",signif(ddx3y.test$estimate,3)), pos = 2)
text(2e5, 6.7e7, paste0("p = ",signif(ddx3y.test$p.value,3)), pos = 2)
dev.off()


#################################
##### C. PRTN3

## DDX3Y
ddx3y.test <- cor.test(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["PRTN3",]), method = "spearman") 

pdf(paste0(Sys.Date(),"_DIA_proteomics_DDX3Y_PRTN3.pdf"))
plot(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["PRTN3",]),
     xlab = "DIA DDX3Y levels (A.U.)",
     ylab = "DIA PRTN3 levels (A.U.)",
     pch = 16,
     col = "tomato",
     xlim = c(0,2e5),
     ylim = c(0,9e7),
     las = 1)
text(2e5, 7.6e7, paste0("Rho = ",signif(ddx3y.test$estimate,3)), pos = 2)
text(2e5, 6.7e7, paste0("p = ",signif(ddx3y.test$p.value,3)), pos = 2)
dev.off()

#################################
##### D. CTSG (CATG in dataset)

ddx3y.test <- cor.test(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["CATG",]), method = "spearman") 

pdf(paste0(Sys.Date(),"_DIA_proteomics_DDX3Y_CTSG.pdf"))
plot(as.numeric(my.mat.prots["DDX3Y",]),as.numeric(my.mat.prots["CATG",]),
     xlab = "DIA DDX3Y levels (A.U.)",
     ylab = "DIA CTSG levels (A.U.)",
     pch = 16,
     col = "tomato",
     xlim = c(0,2e5),
     ylim = c(0,9e7),
     las = 1)
text(2e5, 7.6e7, paste0("Rho = ",signif(ddx3y.test$estimate,3)), pos = 2)
text(2e5, 6.7e7, paste0("p = ",signif(ddx3y.test$p.value,3)), pos = 2)
dev.off()

