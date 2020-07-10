library(pheatmap)
library(RColorBrewer)
library(fields)
library(limma)

#######################################################################################################################################
# https://stackoverflow.com/questions/14271584/r-legend-for-color-density-scatterplot-produced-using-smoothscatter
fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}


normalize_to_BCA <- function(my.matrix, my.vector) {
  my.new.matrix <- matrix(0,nrow(my.matrix),ncol(my.matrix))
  for (i in 1:ncol(my.matrix)) {
    my.new.matrix[,i] <- my.matrix[,i]/my.vector[i]
  }
  return(my.new.matrix)
}
