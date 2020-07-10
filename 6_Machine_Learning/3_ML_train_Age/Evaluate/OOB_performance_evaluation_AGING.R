setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Age_Expression/ML_run/Aging_ML/')
options(stringsAsFactors = F)

source('ML_helper_functions.R')
library(beeswarm)

# 2020-06-02
# Run same framework as sex-dimorphism

########################################################################################
load('../../Features/Feature_engineering/2020-06-01_ALL_ML_features_for_AGING_Expression.RData')

# select columns for ML
my.ALL.features.2group <- unique(my.ALL.features.2group)
rownames(my.ALL.features.2group) <- my.ALL.features.2group$Gene.Name
colnames(my.ALL.features.2group)[1:10]
my.col.exclude <- c("Gene.Name")
my.ALL.features.2group.inc     <- my.ALL.features.2group[,-which(colnames(my.ALL.features.2group) %in% (my.col.exclude))]

# load data partition
load('2020-06-01_Training_Data_AGING.RData')

my.training.data <- my.ALL.features.2group.inc[my.training.idx,]
########################################################################################

#########################################################################
load("2020-06-01_LDA_model_AGING.RData")
load("2020-06-01_NNET_model_AGING.RData")
load("2020-06-01_RF_model_AGING.RData")
load("2020-06-01_SVM_model_AGING.RData")
load("2020-06-02_GBM_model_AGING.RData")
load("2020-06-01_cTree_model_AGING.RData")
load("2020-06-01_Regularized_Logistic_Regression_model_AGING.RData")

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

write.table(my.model.summaries$statistics$balancedAcc, file = paste0(Sys.Date(),"_10_fold_CV_OOB_model_accuracies_summary_AGING.txt"), sep = "\t", quote = F)


my.10cv.data <- my.model.summaries$values
my.median.acc <- apply(my.10cv.data,2,median)

pdf(paste(Sys.Date(),"Machine_Learning_OOB_Accuracy_boxplots_10CV_beeswarm_AGING.pdf", sep =""), height = 5, width = 6)
par(oma=c(0.1,5,0.1,0.1))
boxplot(my.10cv.data[,order(my.median.acc)], las = 1, ylim = c(0.4,1), horizontal = T, xlab = "10 Fold-CV balanced accuracy", outline = F, col = "purple")
beeswarm(my.10cv.data[,order(my.median.acc)], add = T, horizontal = T, pch = 16, cex = 0.5)
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()

#######################
sink(file = paste(Sys.Date(),"Machine_Learning_OOB_performance_evaluation_R_session_Info_AGING.txt", sep =""))
sessionInfo()
sink()


