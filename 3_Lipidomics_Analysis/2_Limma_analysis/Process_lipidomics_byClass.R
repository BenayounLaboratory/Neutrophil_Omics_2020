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

# use lm to determine effect and significance by class

#########################################################
#################    Class analysis     #################
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

####################################
# 4. retrieve annotations
my.lipid.data.vsn <- cbind(rownames(my.lipid.data.vsn),my.lipid.data.vsn)
colnames(my.lipid.data.vsn)[1] <- "Lipid_ID"

my.lipid.data.vsn.annot <- merge(my.lipid.FA.data[,1:6],my.lipid.data.vsn, by = "Lipid_ID")


############################################################################################################
################ by lipid type
my.lipid.data.type <- aggregate(my.lipid.data.vsn.annot[,-c(1:6)], 
                                by = list(my.lipid.data.vsn.annot$Lipid_Family), 
                                FUN = sum)
rownames(my.lipid.data.type) <- my.lipid.data.type$Group.1


pdf(paste(Sys.Date(),"Neutrophils_Aging_Sex_Lipidomics_by_Lipid_Family_BCA_VSN_heatmap_SORTED_BY_CLASS.pdf",sep = "_"), width = 8, height = 5, onefile = F)
pheatmap(my.lipid.data.type[c("PE","PC","LPE","LPC", "SM","CER","DCER","LCER","DAG","TAG","CE","FFA"),-1], 
         scale = 'row', 
         cluster_cols = F,
         cluster_rows = F,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         cellwidth = 20)
dev.off()

# Phospholipids (PE+PC+LPE+LPC)
# Sphingolipids (SM+CER)
# Neutral lipids (DAG+TAG)
# CE Cholesteryl ester

my.lipid.type.list <- list("YF_PE"    = as.numeric(my.lipid.data.type["PE",2:5]  )  ,
                           "YM_PE"    = as.numeric(my.lipid.data.type["PE",10:13])  ,
                           "OF_PE"    = as.numeric(my.lipid.data.type["PE",6:9]  )  ,
                           "OM_PE"    = as.numeric(my.lipid.data.type["PE",14:17])  ,
                           "YF_PC"    = as.numeric(my.lipid.data.type["PC",2:5]  )  ,
                           "YM_PC"    = as.numeric(my.lipid.data.type["PC",10:13])  ,
                           "OF_PC"    = as.numeric(my.lipid.data.type["PC",6:9]  )  ,
                           "OM_PC"    = as.numeric(my.lipid.data.type["PC",14:17])  ,
                           "YF_LPE"   = as.numeric(my.lipid.data.type["LPE",2:5]  )  ,
                           "YM_LPE"   = as.numeric(my.lipid.data.type["LPE",10:13])  ,
                           "OF_LPE"   = as.numeric(my.lipid.data.type["LPE",6:9]  )  ,
                           "OM_LPE"   = as.numeric(my.lipid.data.type["LPE",14:17])  ,
                           "YF_LPC"   = as.numeric(my.lipid.data.type["LPC",2:5]  )  ,
                           "YM_LPC"   = as.numeric(my.lipid.data.type["LPC",10:13])  ,
                           "OF_LPC"   = as.numeric(my.lipid.data.type["LPC",6:9]  )  ,
                           "OM_LPC"   = as.numeric(my.lipid.data.type["LPC",14:17])  ,
                           "YF_SM"    = as.numeric(my.lipid.data.type["SM",2:5]  )  ,
                           "YM_SM"    = as.numeric(my.lipid.data.type["SM",10:13])  ,
                           "OF_SM"    = as.numeric(my.lipid.data.type["SM",6:9]  )  ,
                           "OM_SM"    = as.numeric(my.lipid.data.type["SM",14:17])  ,
                           "YF_CER"   = as.numeric(my.lipid.data.type["CER",2:5]  )  ,
                           "YM_CER"   = as.numeric(my.lipid.data.type["CER",10:13])  ,
                           "OF_CER"   = as.numeric(my.lipid.data.type["CER",6:9]  )  ,
                           "OM_CER"   = as.numeric(my.lipid.data.type["CER",14:17])  ,
                           "YF_LCER"  = as.numeric(my.lipid.data.type["LCER",2:5]  )  ,
                           "YM_LCER"  = as.numeric(my.lipid.data.type["LCER",10:13])  ,
                           "OF_LCER"  = as.numeric(my.lipid.data.type["LCER",6:9]  )  ,
                           "OM_LCER"  = as.numeric(my.lipid.data.type["LCER",14:17])  ,
                           "YF_DCER"  = as.numeric(my.lipid.data.type["DCER",2:5]  )  ,
                           "YM_DCER"  = as.numeric(my.lipid.data.type["DCER",10:13])  ,
                           "OF_DCER"  = as.numeric(my.lipid.data.type["DCER",6:9]  )  ,
                           "OM_DCER"  = as.numeric(my.lipid.data.type["DCER",14:17])  ,
                           "YF_DAG"   = as.numeric(my.lipid.data.type["DAG",2:5]  )  ,
                           "YM_DAG"   = as.numeric(my.lipid.data.type["DAG",10:13])  ,
                           "OF_DAG"   = as.numeric(my.lipid.data.type["DAG",6:9]  )  ,
                           "OM_DAG"   = as.numeric(my.lipid.data.type["DAG",14:17])  ,
                           "YF_TAG"   = as.numeric(my.lipid.data.type["TAG",2:5]  )  ,
                           "YM_TAG"   = as.numeric(my.lipid.data.type["TAG",10:13])  ,
                           "OF_TAG"   = as.numeric(my.lipid.data.type["TAG",6:9]  )  ,
                           "OM_TAG"   = as.numeric(my.lipid.data.type["TAG",14:17])  ,
                           "YF_CE"    = as.numeric(my.lipid.data.type["CE",2:5]  )  ,
                           "YM_CE"    = as.numeric(my.lipid.data.type["CE",10:13])  ,
                           "OF_CE"    = as.numeric(my.lipid.data.type["CE",6:9]  )  ,
                           "OM_CE"    = as.numeric(my.lipid.data.type["CE",14:17])  ,
                           "YF_FFA"   = as.numeric(my.lipid.data.type["FFA",2:5]  )  ,
                           "YM_FFA"   = as.numeric(my.lipid.data.type["FFA",10:13])  ,
                           "OF_FFA"   = as.numeric(my.lipid.data.type["FFA",6:9]  )  ,
                           "OM_FFA"   = as.numeric(my.lipid.data.type["FFA",14:17]) 
) 

pdf(paste(Sys.Date(),"Neutrophils_Aging_Sex_Lipidomics_Lipid_types_Boxplot.pdf",sep = "_"), width = 15, height = 10)
boxplot(my.lipid.type.list,
        col = c("deeppink","deepskyblue","deeppink4","deepskyblue4"),
        las = 2, 
        ylab = "Amount (A.U.)", 
        #ylim = c(0,1000),
        main = "Lipids by class",
        log = 'y')
dev.off()


###### statistical test
my.class.lips <- list("PE"   =   1:4,
                      "PC"   =   5:8,
                      "LPE"  =  9:12,
                      "LPC"  = 13:16,
                      "SM"   = 17:20,
                      "CER"  = 21:24,
                      "LCER" = 25:28,
                      "DCER" = 29:32,
                      "DAG"  = 33:36,
                      "TAG"  = 37:40,
                      "CE"   = 41:44,
                      "FFA"  = 45:48)


my.lm.sigs <- data.frame("Sex_Effect_Coefficient" = rep(NA,length(my.class.lips)),
                         "Age_Effect_Coefficient" = rep(NA,length(my.class.lips)),
                         "Sex_Effect_pval"        = rep(NA,length(my.class.lips)) ,
                         "Age_Effect_pval"        = rep(NA,length(my.class.lips)) )
rownames(my.lm.sigs) <- names(my.class.lips)

for (i in 1:length(my.class.lips)) {
        my.cur.lip.data <- data.frame("Lip"  = unlist(my.lipid.type.list[my.class.lips[[i]]]),
                                      "Sex" = c(rep("F",4),rep("M",4),rep("F",4),rep("M",4)),
                                      "Age" = c(rep(4,8),rep(20,8)))
        my.res <- summary(lm(Lip ~ Sex + Age, data = my.cur.lip.data))
        
        my.lm.sigs[i,1:2] <- as.numeric(my.res$coefficients[-1,1]) # coefficient estimate
        my.lm.sigs[i,3:4] <- as.numeric(my.res$coefficients[-1,4]) # coefficient p-value
        
}

my.lm.sigs$Sex_Effect_fdr <- p.adjust(my.lm.sigs$Sex_Effect_p, method = "BH")
my.lm.sigs$Age_Effect_fdr <- p.adjust(my.lm.sigs$Age_Effect_p, method = "BH")

write.table(my.lm.sigs,file = paste0(Sys.Date(), "_Lipid_by_Class_LM_Significance.txt"), sep = "\t", quote = F)

#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()



