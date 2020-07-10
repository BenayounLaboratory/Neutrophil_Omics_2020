setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Age_Expression/ML_run/Aging_ML/')
options(stringsAsFactors = F)

source('ML_helper_functions.R')

load('../../Features/Feature_engineering/2020-06-01_ALL_ML_features_for_AGING_Expression.RData')

# Run same framework as sex-dimorphism

########################################################################################
# select columns for ML
my.ALL.features.2group <- unique(my.ALL.features.2group)

rownames(my.ALL.features.2group) <- my.ALL.features.2group$Gene.Name
colnames(my.ALL.features.2group)[1:10]
my.col.exclude <- c("Gene.Name")

my.ALL.features.2group.inc     <- my.ALL.features.2group[,-which(colnames(my.ALL.features.2group) %in% (my.col.exclude))]

########################################################################################
# get data partition
set.seed(123456789)
my.training.idx <- createDataPartition(my.ALL.features.2group.inc$AGING, p=0.67, list=FALSE) # 2/3 for training

my.training.data <- my.ALL.features.2group.inc[my.training.idx,]

save(my.training.data, my.training.idx, file = paste0(Sys.Date(),"_Training_Data_AGING.RData"))

########################################################################################
##### 1. run SVM model

# use 10-fold cross-validation to build the model
my.ctrl.opt.svm <- trainControl(method = "cv",
                                number = 10,
                                allowParallel=TRUE,
                                verbose=F,
                                summaryFunction=getTwoClassBalancedAccuracy,
                                classProbs = TRUE)

# train model with caret train function
my.svm.fit       <- train( AGING ~ .,
                           data = my.training.data,
                           method="svmRadial",
                           trControl = my.ctrl.opt.svm,
                           tuneLength = 10,
                           metric="balancedAcc")

save(my.svm.fit, file = paste0(Sys.Date(),"_SVM_model_AGING.RData"))


########################################################################################
##### 2. run Sparse LDA model

# use 10-fold cross-validation to build the model
my.ctrl.opt.lda <- trainControl(method = "cv",
                                number = 10,
                                allowParallel=TRUE,
                                verbose=F,
                                summaryFunction=getTwoClassBalancedAccuracy,
                                classProbs = TRUE)

fineGrid.LDA <-  expand.grid(NumVars = seq(15,200,15), 
                             lambda = c(0, 0.01, 0.1, 1, 10, 100))

# train model with caret train function
my.lda.fit       <- train( AGING ~ .,
                           data = my.training.data,
                           method="sparseLDA",
                           trControl = my.ctrl.opt.lda,
                           tuneGrid = fineGrid.LDA,
                           metric="balancedAcc")

save(my.lda.fit, file = paste0(Sys.Date(),"_LDA_model_AGING.RData"))

########################################################################################
##### 3. run NNET model
set.seed(123456789)

# use 10-fold cross-validation to build the model
my.ctrl.opt.nnet <- trainControl(method = "cv",
                                 number = 10,
                                 allowParallel=TRUE,
                                 verbose=F,
                                 summaryFunction=getTwoClassBalancedAccuracy,
                                 classProbs = TRUE)

# train model with caret train function
my.nnet.fit <- train( AGING ~ .,
                      data = my.training.data,
                      method="pcaNNet",
                      trControl = my.ctrl.opt.nnet,
                      tuneLength = 10,
                      metric="balancedAcc")

save(my.nnet.fit, file = paste0(Sys.Date(),"_NNET_model_AGING.RData"))


########################################################################################
##### 4. run Regularized Logistic Regression model
set.seed(123456789)

# use 10-fold cross-validation to build the model
my.ctrl.opt.logreg <- trainControl(method = "cv",
                                   number = 10,
                                   allowParallel=TRUE,
                                   verbose=F,
                                   summaryFunction=getTwoClassBalancedAccuracy,
                                   classProbs = TRUE)

# train model with caret train function
my.logreg.fit <- train( AGING ~ .,
                        data = my.training.data,
                        method="regLogistic",
                        trControl = my.ctrl.opt.logreg,
                        tuneLength = 10,
                        metric="balancedAcc")

save(my.logreg.fit, file = paste0(Sys.Date(),"_Regularized_Logistic_Regression_model_AGING.RData"))


########################################################################################
##### 5. run RF model
set.seed(123456789)

# use 10-fold cross-validation to build the model
my.ctrl.opt.rf <- trainControl(method = "cv",
                               number = 10,
                               allowParallel=TRUE,
                               verbose=F,
                               summaryFunction=getTwoClassBalancedAccuracy,
                               classProbs = TRUE)

fineGrid.rf <- expand.grid(mtry = seq(3,21,3))

# train model with caret train function
my.rf.fit <- train( AGING ~ .,
                    data = my.training.data,
                    method="rf",
                    importance=TRUE,
                    trControl = my.ctrl.opt.rf,
                    tuneGrid = fineGrid.rf,
                    metric="balancedAcc"
)

save(my.rf.fit, file = paste0(Sys.Date(),"_RF_model_AGING.RData"))

########################################################################################
##### 6. run Conditional Inference Tree model
set.seed(123456789)

# use 10-fold cross-validation to build the model
my.ctrl.opt.ctree <- trainControl(method = "cv",
                                  number = 10,
                                  allowParallel=TRUE,
                                  verbose=F,
                                  summaryFunction=getTwoClassBalancedAccuracy,
                                  classProbs = TRUE)

# train model with caret train ufnction
my.ctree.fit       <- train( AGING ~ .,
                             data = my.training.data,
                             method="ctree",
                             trControl = my.ctrl.opt.ctree,
                             tuneLength = 10,
                             metric="balancedAcc")

save(my.ctree.fit, file = paste0(Sys.Date(),"_cTree_model_AGING.RData"))


########################################################################################
##### 7. run GBM model
set.seed(123456789)

# use 10-fold cross-validation to build the model
my.ctrl.opt.gbm <- trainControl(method = "cv",
                                number = 10,
                                allowParallel=TRUE,
                                verbose=F,
                                summaryFunction=getTwoClassBalancedAccuracy,
                                classProbs = TRUE)

fineGrid.gbm <- expand.grid(n.trees=seq(6000,14000,4000), 
                            interaction.depth = seq(6,12,3),
                            shrinkage = 0.001, 
                            n.minobsinnode = 20)


# train model with caret train function
my.gbm.fit <- train( AGING ~ .,
                     data = my.training.data,
                     method="gbm",
                     trControl = my.ctrl.opt.gbm,
                     tuneGrid = fineGrid.gbm,
                     metric="balancedAcc")


save(my.gbm.fit, file = paste0(Sys.Date(),"_GBM_model_AGING.RData"))

#######################
sink(file = paste(Sys.Date(),"Machine_Learning_R_session_Info_AGING.txt", sep =""))
sessionInfo()
sink()

