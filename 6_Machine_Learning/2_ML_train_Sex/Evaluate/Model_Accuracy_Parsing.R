setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/Machine_Learning/ML_Predict_Sex_Dim_Expression/ML_run/Testing_accuracy/')
options(stringsAsFactors = F)

source('accuracy_parsing_functions.R')

# Get testing accuracy on trained models
# calculate variable importance TF

########################################################################################
# Load and reformart trainign data for ML
load('../../Features/Feature_engineering/2020-05-20_ALL_ML_features_for_SexDim_Expression.RData')
load('../Sex_Dim_ML/2020-05-20_Training_Data.RData')

rownames(my.ALL.features.2group) <- my.ALL.features.2group$Gene.Name
colnames(my.ALL.features.2group)[1:10]
my.col.exclude <- c("Gene.Name","FDR","logFC")

my.ALL.features.2group.inc     <- my.ALL.features.2group[,-which(colnames(my.ALL.features.2group) %in% (my.col.exclude))]

my.testing.data <- my.ALL.features.2group.inc[-my.training.idx,]

########################################################################################
# Load models
load("../Sex_Dim_ML/2020-05-20_LDA_model.RData")
load("../Sex_Dim_ML/2020-05-20_NNET_model.RData")
load("../Sex_Dim_ML/2020-05-20_RF_model.RData")
load("../Sex_Dim_ML/2020-05-20_SVM_model.RData")
load("../Sex_Dim_ML/2020-05-20_GBM_model.RData")
load("../Sex_Dim_ML/2020-05-20_cTree_model.RData")
load("../Sex_Dim_ML/2020-05-20_Regularized_Logistic_Regression_model.RData")

Acc.svm.fit     <- get_acc_metrics("Neutrophil_Sex_Dim_SVM"                  , my.svm.fit      , my.testing.data)
Acc.lda.fit     <- get_acc_metrics("Neutrophil_Sex_Dim_LDA"                  , my.lda.fit      , my.testing.data)
Acc.nnet.fit    <- get_acc_metrics("Neutrophil_Sex_Dim_NNET"                 , my.nnet.fit     , my.testing.data)
Acc.rf.fit      <- get_acc_metrics("Neutrophil_Sex_Dim_RF"                   , my.rf.fit       , my.testing.data)
Acc.gbm.fit     <- get_acc_metrics("Neutrophil_Sex_Dim_GBM"                  , my.gbm.fit      , my.testing.data)
Acc.ctree.fit   <- get_acc_metrics("Neutrophil_Sex_Dim_cTree"                , my.ctree.fit    , my.testing.data)
Acc.logreg.fit <- get_acc_metrics("Neutrophil_Sex_Dim_Regularized_Logistic"  , my.logreg.fit   , my.testing.data)

my.testing.accuracies       <- rbind(Acc.svm.fit     ,
                                     Acc.lda.fit     ,
                                     Acc.nnet.fit    ,
                                     Acc.rf.fit      ,
                                     Acc.gbm.fit     ,
                                     Acc.ctree.fit   ,
                                     Acc.logreg.fit   )

rownames(my.testing.accuracies) <- c("SVM" ,
                                     "LDA" ,
                                     "NNET" ,
                                     "RF" ,
                                     "GBM" ,
                                     "cTree" ,
                                     "LogReg" )

colnames(my.testing.accuracies) <- c("Balanced_accuracy","AUC")

my.testing.accuracies.srt <- my.testing.accuracies[c("RF"    ,
                                                     "GBM"   ,
                                                     "LogReg",
                                                     "NNET"  ,
                                                     "SVM"   ,
                                                     "LDA"   ,
                                                     "cTree" ) ,]
my.testing.accuracies.srt <- data.frame(my.testing.accuracies.srt)

pdf(paste0(Sys.Date(),"_Testing_Balanced_accuracy_dotplot.pdf"))
dotchart(rev(my.testing.accuracies.srt$Balanced_accuracy),
         labels = rev(rownames (my.testing.accuracies.srt)),
         las = 1,
         xlim = c(0.4,1),
         col= "purple",
         pt.cex = 3,
         pch = 16,
         main = "Balanced Accuracy")
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()

pdf(paste0(Sys.Date(),"_Testing_AUC_dotplot.pdf"))
dotchart(rev(my.testing.accuracies.srt$AUC),
         labels = rev(rownames (my.testing.accuracies.srt)),
         las = 1,
         xlim = c(0.4,1),
         col= "salmon1",
         pt.cex = 3,
         pch = 16,
         main = "AUC")
abline(v = 0.5, col = "red", lty = "dashed")
dev.off()

#### get heatmap
library('pheatmap')

my.accuracy.colors <- c("floralwhite","lightsalmon","indianred1","firebrick3","firebrick","firebrick4")

my.testing.accuracies.srt.heat <- rbind(my.testing.accuracies.srt,
                                        c(0.5,0.5)               ,
                                        c(1,1)                     )
rownames(my.testing.accuracies.srt.heat)[8:9] <- c("Random","Perfect")

pdf(paste0(Sys.Date(),"_Testing_Accuracy_AUC_heatmap.pdf"))
pheatmap(my.testing.accuracies.srt.heat,
         cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(my.accuracy.colors)(50),
         cellwidth = 15,
         cellheight = 15)
dev.off()


########################################################################################
# Get important features
#### https://topepo.github.io/caret/variable-importance.html

# feature importance on tree-based models 
my.rf.Imp    <- varImp(my.rf.fit    , useModel = TRUE , scale = TRUE)$importance
my.gbm.Imp   <- varImp(my.gbm.fit   , useModel = TRUE , scale = TRUE)$importance

my.var.Imps <- data.frame(matrix(0,nrow(my.rf.Imp),5))
colnames(my.var.Imps) <- c("Feature","VarImp_GBM","VarImp_RF","Rank_GBM","Rank_RF")
my.var.Imps$Feature <- rownames(my.rf.Imp)

# parse varImp calculation
for ( i in 1:nrow(my.rf.Imp)) {
        
        my.rf.ix  <- i
        my.gbm.ix <- rownames(my.gbm.Imp) %in% rownames(my.rf.Imp)[i]

        my.var.Imps$VarImp_RF[i]    <-  apply(my.rf.Imp,1,mean)[my.rf.ix]
        my.var.Imps$VarImp_GBM[i]   <-  my.gbm.Imp$Overall[my.gbm.ix]

}

# now get ranks by decreasing order
my.sort.rf   <- sort(my.var.Imps$VarImp_RF  , index.return = T, decreasing = TRUE)
my.sort.gbm  <- sort(my.var.Imps$VarImp_GBM , index.return = T, decreasing = TRUE)

my.var.Imps$Rank_RF[my.sort.rf$ix]    <- 1:nrow(my.var.Imps)
my.var.Imps$Rank_GBM[my.sort.gbm$ix]  <- 1:nrow(my.var.Imps)

my.var.Imps$RANK_Product               <-  my.var.Imps$Rank_RF * my.var.Imps$Rank_GBM

save(my.var.Imps, file = paste0(Sys.Date(),"_Variable_Importance_Parsing_RF_GBM_CARET.RData"))
write.table(my.var.Imps,file = paste0(Sys.Date(),"_Variable_Importance_Parsing_RF_GBM_CARET.txt"), quote = F, row.names = F, sep = "\t")

# feature importance on tree based models
my.rf.varimps.native <- my.rf.fit$finalModel$importance[,3:4]
save(my.rf.varimps.native, file = paste0(Sys.Date(),"_Variable_Importance_RF_native.RData"))
########################################################################################


########################################################################################
# plot combined accuracy
library(ggplot2)

my.col.palette.acc <- colorRampPalette(rev(c("#FF9999","indianred1","firebrick1","firebrick3","firebrick4")))(1000)

my.rank.sort <- sort (my.var.Imps$RANK_Product, decreasing = F, index.return = T)
my.var.Imps.sorted <- my.var.Imps[rev(my.rank.sort$ix[1:20]),]

# format for ggplot; initialize with RF
my.rf.results   <- data.frame('Feature_name' =  my.var.Imps.sorted$Feature    ,
                              'VarImp'       =  my.var.Imps.sorted$VarImp_RF ,
                              'Rank'         =  my.var.Imps.sorted$Rank_RF   ,
                              'Model'        = "RF"
)

my.gbm.results <- data.frame('Feature_name'  =  my.var.Imps.sorted$Feature   ,
                             'VarImp'       =  my.var.Imps.sorted$VarImp_GBM ,
                             'Rank'         =  my.var.Imps.sorted$Rank_GBM   ,
                             'Model'        = "GBM"
)

my.complete.results <- rbind(my.rf.results,
                             my.gbm.results)

# will make Importance the size and Rank the color

# to preserve the wanted order
my.complete.results$Feature_name <- factor(my.complete.results$Feature_name, levels = unique(my.complete.results$Feature_name))

my.pdfname <-paste(Sys.Date(),"RF_GBM_Models_VarImps_balloon_plot_top_20.pdf", sep="_")

pdf(my.pdfname, onefile=F, height = 9, width=5)
my.plot <- ggplot(my.complete.results,aes(x=Model,y=Feature_name,colour=Rank,size=VarImp))+ theme_bw()+ geom_point(shape = 16) 
my.plot <- my.plot + ggtitle("Model Feature Importance") + labs(x = "Tissue", y = "Feature")
my.plot <- my.plot + scale_colour_gradientn(colours = my.col.palette.acc,space = "Lab", na.value = "grey50", guide = "colourbar")
print(my.plot)
dev.off()  






