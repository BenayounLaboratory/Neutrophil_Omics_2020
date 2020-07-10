# Data processing BB Mouse Neutrophils
# Imputation
tmp <- lipid_C_Ber_data
tmp[is.na(tmp)] <- 1
tmp <- log2(tmp)
tmp[tmp == 0] <- NA
tmp <- tmp[rowSums(!is.na(tmp)) > 0.66*ncol(tmp),]
lipid_class <- c("CE.[[:digit:]]","^CER","DAG",".CER","FFA","LPC","LPE","^PC","^PE","SM","TAG")

tmp2 <- data.frame()
pdf("NA_imputation_lipid_C_Shay_data_1.8sd.pdf")
for (j in 1:length(lipid_class)){
  sub <- tmp[grep(lipid_class[j], rownames(tmp)),]  
  for (i in 1:ncol(tmp)){
    y <- which(is.na(sub[,i]) == T)
    sub[which(is.na(sub[,i]) == T),i] <- rnorm(n = length(which(is.na(sub[,i]) == T)), 
                                               mean = mean(as.numeric(sub[,i]),na.rm = T)-(1.8*sd(sub[,i],na.rm = T)), 
                                               sd = sd(sub[,i],na.rm = T)*0.3)
    ggdata <- data.frame(Log2_Int = sub[,i])
    ggdata$Condition <- '1'
    ggdata$Condition[y] <- '0'
    p <- ggplot(ggdata, aes(x = Log2_Int, fill = Condition)) + geom_bar(color="black",binwidth = 0.5) + 
      theme_bw() + ggtitle(colnames(tmp)[i]) + scale_fill_manual(values = c("orange", "grey"))
    print(p)
  }
  tmp2 <- rbind(tmp2,sub)
}
dev.off()
tmp2[is.na(tmp2)]

tmp3 <- 2^tmp2
Lipidomics <- tmp3
rownames(Lipidomics) <- rownames(tmp)
Lipidomics <- Lipidomics[,order(colnames(Lipidomics), decreasing = F)]
Lipidomics <- data.frame(Lipidomics)
write.csv(Lipidomics, "Lipidomics_lipid_C_Ber_data.csv")