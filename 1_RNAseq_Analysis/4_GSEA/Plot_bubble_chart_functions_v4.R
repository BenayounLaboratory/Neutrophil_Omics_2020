library(ggplot2) 
library(scales) 
theme_set(theme_bw())


get_sex_bubble_plot <- function(my.geneset.name, my.gsea.file, max.path.plot = 20){
  
  my.gsea.sex <- read.csv(my.gsea.file, sep = "\t", header = T)
  
  my.gsea.pos <- my.gsea.sex[my.gsea.sex$nes > 0,]
  my.gsea.neg <- my.gsea.sex[my.gsea.sex$nes < 0,]
  
  my.pos.sort <- sort(my.gsea.pos$nes, index.return = T, decreasing = T) # largest value is top (positive)
  my.neg.sort <- sort(my.gsea.neg$nes, index.return = T, decreasing = F) # largest value is top (negative)
  
  if ( (nrow(my.gsea.pos) > round(max.path.plot/2)) && (nrow(my.gsea.neg) > round(max.path.plot/2)) ) {
    
    # if enough on both sides
    my.gsea.sex.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:round(max.path.plot/2)],],
                           my.gsea.neg[my.neg.sort$ix[1:round(max.path.plot/2)],])
    
  } else {
    
    # if not enough on both sides
    my.gsea.sex.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:(min(round(max.path.plot/2),nrow(my.gsea.pos)))],],
                           my.gsea.neg[my.neg.sort$ix[1:(min(round(max.path.plot/2),nrow(my.gsea.neg)))],])
    
  }
  
  # create -log10 FDR for plotting
  my.gsea.sex.2$minlog10fdr  <- -log10(my.gsea.sex.2$fdr + 1e-30)
  my.gsea.sex.2$Description  <- rownames(my.gsea.sex.2)
  
  my.sorting <- sort(my.gsea.sex.2$minlog10fdr, index.return = T, decreasing = T)
  my.gsea.sex.sorted <- my.gsea.sex.2[my.sorting$ix,]
  
  # create and preserve wanted display order
  my.gsea.sex.sorted$sex <- ifelse(my.gsea.sex.sorted$nes < 0, "Male", "Female")  # male/female avg flag
  my.gsea.sex.sorted <- my.gsea.sex.sorted[order(my.gsea.sex.sorted$sex), ]
  
  my.max.char <- max(nchar(my.gsea.sex.sorted$Description))
  
  my.gsea.sex.sorted$Description <- factor(my.gsea.sex.sorted$Description, levels = rev(unique(my.gsea.sex.sorted$Description)))
  my.gsea.sex.sorted$CellType <- factor(rep("Neutrophils",length(  my.gsea.sex.sorted$Description)))
  
  # Female/Male color scale
  my.max <- max(my.gsea.sex.sorted$nes)
  my.min <- min(my.gsea.sex.sorted$nes)
  
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector.sex <- c("deepskyblue","lightskyblue","lightskyblue1","lightcyan","white","lavenderblush","plum1","orchid1","deeppink")
  
  my.plot <- ggplot(my.gsea.sex.sorted,aes(x=CellType,y=Description,colour=nes,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle("gsea Analysis") + labs(x = "-log10(pvalue)", y = "")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.sex, na.value = "grey50", guide = "colourbar", values = my.scaled)
  my.plot
  
  my.pdfname <- paste(Sys.Date(),"GSEA_sex_BALLOON_plot",my.geneset.name,"top", nrow(my.gsea.sex.sorted),"significant_pathways.pdf", sep="_")
  
  pdf(my.pdfname, onefile=T, height = 6, width=max(4,my.max.char/8) )
  print(my.plot)
  dev.off()  
  
  
}



##########################################################################################
get_age_bubble_plot <- function(my.geneset.name, my.gsea.file, max.path.plot = 20){
  
  my.gsea.age <- read.csv(my.gsea.file, sep = "\t", header = T)
  
  my.gsea.pos <- my.gsea.age[my.gsea.age$nes > 0,]
  my.gsea.neg <- my.gsea.age[my.gsea.age$nes < 0,]
  
  my.pos.sort <- sort(my.gsea.pos$nes, index.return = T, decreasing = T) # largest value is top (positive)
  my.neg.sort <- sort(my.gsea.neg$nes, index.return = T, decreasing = F) # largest value is top (negative)
  
  if ( (nrow(my.gsea.pos) > round(max.path.plot/2)) && (nrow(my.gsea.neg) > round(max.path.plot/2)) ) {
    
    # if enough on both sides
    my.gsea.age.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:round(max.path.plot/2)],],
                           my.gsea.neg[my.neg.sort$ix[1:round(max.path.plot/2)],])
    
  } else {
    
    # if not enough on both sides
    my.gsea.age.2 <- rbind(my.gsea.pos[my.pos.sort$ix[1:(min(round(max.path.plot/2),nrow(my.gsea.pos)))],],
                           my.gsea.neg[my.neg.sort$ix[1:(min(round(max.path.plot/2),nrow(my.gsea.neg)))],])
    
  }
  
  # create -log10 FDR for plotting
  my.gsea.age.2$minlog10fdr  <- -log10(my.gsea.age.2$fdr + 1e-30)
  my.gsea.age.2$Description  <- rownames(my.gsea.age.2)
  
  my.sorting <- sort(my.gsea.age.2$minlog10fdr, index.return = T, decreasing = T)
  my.gsea.age.sorted <- my.gsea.age.2[my.sorting$ix,]
  
  # create and preserve wanted display order
  my.gsea.age.sorted$direction <- ifelse(my.gsea.age.sorted$nes > 0, "Up", "Down")  # Up down avg flag
  my.gsea.age.sorted <- my.gsea.age.sorted[rev(order(my.gsea.age.sorted$direction)), ]
  
  my.max.char <- max(nchar(my.gsea.age.sorted$Description))
  
  my.gsea.age.sorted$Description <- factor(my.gsea.age.sorted$Description, levels = rev(unique(my.gsea.age.sorted$Description)))
  my.gsea.age.sorted$CellType <- factor(rep("Neutrophils",length(  my.gsea.age.sorted$Description)))
  
  # Female/Male color scale
  my.max <- max(my.gsea.age.sorted$nes)
  my.min <- min(my.gsea.age.sorted$nes)
  
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector.age <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
  
  my.plot <- ggplot(my.gsea.age.sorted,aes(x=CellType,y=Description,colour=nes,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
  my.plot <- my.plot + ggtitle("gsea Analysis") + labs(x = "-log10(pvalue)", y = "")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.age, na.value = "grey50", guide = "colourbar", values = my.scaled)
  my.plot
  
  my.pdfname <- paste(Sys.Date(),"GSEA_aging_BALLOON_plot",my.geneset.name,"top", nrow(my.gsea.age.sorted),"significant_pathways.pdf", sep="_")
  
  pdf(my.pdfname, onefile=T, height = 6, width=max(4,my.max.char/8))
  print(my.plot)
  dev.off()  
  
  
}
