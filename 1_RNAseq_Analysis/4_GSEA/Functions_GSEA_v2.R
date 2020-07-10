library(DESeq2)
library(phenoTest)
library(qusage)   

#################################################################
run_enrich <- function(my.matrix, data.name, my.fdr = 0.05, my.ontology, my.ontology.name) {
  
  # set seed to stabilize output (add 2020/5/20)
  set.seed(123456789)
  
  # run phenotest GSEA
  gsea.data <- gsea( x         =  my.matrix    , 
                     gsets     =  my.ontology , 
                     mc.cores  =  2            , 
                     logScale  =  FALSE        , 
                     B         =  10000         , 
                     minGenes  =  5            , 
                     maxGenes  =  5000         )
  
  my.summary <- data.frame(summary(gsea.data))
  my.sig.path.num <- sum(my.summary$fdr < my.fdr )
  
  # write results to file
  my.outfile <- paste(Sys.Date(),data.name, my.ontology.name, "FDR", 100*my.fdr, "Phenotest_GSEA_Analysis_table", my.sig.path.num,"significant.txt", sep = "_")
  write.table(my.summary[my.summary$fdr < my.fdr,], file = my.outfile, quote = F, sep = "\t")
  
  return(my.summary)
}

