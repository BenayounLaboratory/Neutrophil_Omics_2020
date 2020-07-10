library(pROC)
library(caret)
library(kernlab)         # Radial SVM
library(nnet)            # PCA Neural Net
library(randomForest)    # random forest
library(gbm)             # GBM
library(sparseLDA)       # penalized LDA
library(e1071)
library(party)
library(caTools)
library(LiblineaR)
#library(doMC)

## DEBUG
# my.model.name <- "Neutrophil_Aging_SVM"
# my.mod.fit    <- my.svm.fit
# my.testing    <- my.testing.data


get_acc_metrics <- function(my.model.name, my.mod.fit, my.testing) {
  
  my.mod.preds  <- predict(my.mod.fit, my.testing)
  my.confus.mat <- confusionMatrix(my.mod.preds,my.testing$AGING) # prediction, then ref
  
  my.bal.acc <- my.confus.mat$byClass["Balanced Accuracy"]
  my.sens    <- my.confus.mat$byClass["Sensitivity"]
  my.spe     <- my.confus.mat$byClass["Specificity"]
  my.prec    <- my.confus.mat$byClass["Precision"]
  
  # get ROC for 2-group prediction
  my.roc <- get_roc(my.mod.fit, my.testing)
  my.auc <- auc(my.roc)
  
  my.filename <- paste(Sys.Date(),my.model.name,"ROC_testing_classification_AGING.pdf",sep="_")
  
  pdf(my.filename)
  plot(1 - my.roc$specificities, my.roc$sensitivities,type='l',col="firebrick", xlab = "1-specificity (False positive rate)", ylab = "sensitivity (True positive rate)")
  abline(0,1,col="grey",lty='dashed')
  text(0.8,0.2,paste("AUC = ",round(my.auc,3),sep=""))
  dev.off()
  
  
  write.table(data.frame(balacc      = my.bal.acc, 
                         AUC         = my.auc,
                         specificity = my.spe,
                         sensitivity = my.sens,
                         precision   = my.prec), 
              row.names = F,
              file= paste(Sys.Date(),my.model.name,"testing_metrics_classification_AGING.txt", sep="_"))
	
	return (c(my.bal.acc,my.auc))
}




#######################################################################################################################################
# get ROC predictions
get_roc <- function(my.mod.fit,my.testing) {
  my.testing.preds.probs <- predict(my.mod.fit,my.testing,type="prob")

  my.roc <- roc(my.testing$AGING, my.testing.preds.probs$DOWN) # Calculate ROC curve.
  
  return(my.roc)
  
}
