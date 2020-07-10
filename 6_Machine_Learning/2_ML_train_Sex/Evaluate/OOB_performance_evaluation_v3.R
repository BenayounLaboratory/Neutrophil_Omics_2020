setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/ML_run/Sex_Dim_ML')
options(stringsAsFactors = F)

source('ML_helper_functions.R')
library(beeswarm)

# Evaluate OOB/CV performance


########################################################################################
load('../../Features/Feature_engineering/2020-05-20_ALL_ML_features_for_SexDim_Expression.RData')

# select columns for ML
my.ALL.features.2group <- unique(my.ALL.features.2group)
rownames(my.ALL.features.2group) <- my.ALL.features.2group$Gene.Name
colnames(my.ALL.features.2group)[1:10]
my.col.exclude <- c("Gene.Name")
my.ALL.features.2group.inc     <- my.ALL.features.2group[,-which(colnames(my.ALL.features.2group) %in% (my.col.exclude))]

# load data partition
load('2020-05-20_Training_Data.RData')

my.training.data <- my.ALL.features.2group.inc[my.training.idx,]
########################################################################################

#########################################################################
load("2020-05-20_LDA_model.RData")
load("2020-05-20_NNET_model.RData")
load("2020-05-20_RF_model.RData")
load("2020-05-20_SVM_model.RData")
load("2020-05-20_GBM_model.RData")
load("2020-05-20_cTree_model.RData")
load("2020-05-20_Regularized_Logistic_Regression_model.RData")

results <- resamples(list("SVM"        = my.svm.fit    ,
                          "LDA"        = my.lda.fit    ,
                          "NNET"       = my.nnet.fit   ,
                          "RF"         = my.rf.fit     ,
                          "GBM"        = my.gbm.fit    ,
                          "cTree"      = my.ctree.fit  ,
                          "LogReg"     = my.logreg.fit
                          ))

# summary of model differences
my.model.summaries <- summary(results)

write.table(my.model.summaries$statistics$balancedAcc, file = paste0(Sys.Date(),"_10_fold_CV_OOB_model_accuracies_summary.txt"), sep = "\t", quote = F)

my.10cv.data <- my.model.summaries$values
my.median.acc <- apply(my.10cv.data,2,median)

pdf(paste(Sys.Date(),"Machine_Learning_OOB_Accuracy_boxplots_10CV_beeswarm.pdf", sep =""), height = 5, width = 6)
par(oma=c(0.1,5,0.1,0.1))
boxplot(my.10cv.data[,order(my.median.acc)], las = 1, ylim = c(0.4,1), horizontal = T, xlab = "10 Fold-CV balanced accuracy", outline = F, col = "purple")
beeswarm(my.10cv.data[,order(my.median.acc)], add = T, horizontal = T, pch = 16, cex = 0.5)
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()

#######################
sink(file = paste(Sys.Date(),"Machine_Learning_OOB_performance_evaluation_R_session_Info.txt", sep =""))
sessionInfo()
sink()


