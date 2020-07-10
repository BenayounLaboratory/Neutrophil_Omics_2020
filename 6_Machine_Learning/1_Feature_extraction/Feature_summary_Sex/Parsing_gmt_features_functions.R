library('GSEABase')
library(clusterProfiler)

######### Process GMT files for use with ClsuterProfiler
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


