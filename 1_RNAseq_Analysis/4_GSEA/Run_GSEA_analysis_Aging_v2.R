setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/GSEA')
options(stringsAsFactors = F)

source('Functions_GSEA_v2.R')

load('2020-05-20_GeneSetCollections_for_Phenotest_GSEA.RData') 
load('2020-06-04Expressed_TF_summary_GeneSetCollections_for_Phenotest_GSEA.RData')

######################## A. Load DEseq2 results for analysis ######################## 
# load DEseq2 results
load('../DEseq2/STAR/2020-05-21_Neutrophils_Aging_BOTH.RData')

my.ntph.age           <- data.frame(res.age)
my.ntph.age$GeneName  <- rownames(my.ntph.age )


######################## B. Prepare GeneLists using DEseq2 t-statistic to rank genes ######################## 
ntph.age.geneList         <- my.ntph.age$stat
names(ntph.age.geneList)  <- my.ntph.age$GeneName
ntph.age.geneList         <- sort(ntph.age.geneList , decreasing = TRUE)


######################## C. Gene Set Enrichment Analysis ######################## 
ntph.age.c2.cp.kegg         <-  run_enrich(ntph.age.geneList, "Neutrophil_Aging_Effect",  my.fdr = 0.05,  my.ontology =  Sym.c2.cp.kegg            , my.ontology.name = "MSigDB_KEGG"               )  
ntph.age.c2.cp.reactome     <-  run_enrich(ntph.age.geneList, "Neutrophil_Aging_Effect",  my.fdr = 0.05,  my.ontology =  Sym.c2.cp.reactome        , my.ontology.name = "MSigDB_Reactome"           )  
ntph.age.ENS.GO.ALL         <-  run_enrich(ntph.age.geneList, "Neutrophil_Aging_Effect",  my.fdr = 0.05,  my.ontology =  Sym.ENS.GO.ALL            , my.ontology.name = "ENS.GO.ALL"                )  
# ntph.age.TF.summary         <-  run_enrich(ntph.age.geneList, "Neutrophil_Aging_Effect",  my.fdr = 0.05,  my.ontology =  Sym.TF.target.summary     , my.ontology.name = "TF_Targets_ALL"                )  

ntph.sex.TF.summary.exp     <-  run_enrich(ntph.age.geneList, "Neutrophil_Aging_Effect",  my.fdr = 0.05,  my.ontology =  Sym.TF.target.summary.exp , my.ontology.name = "TF_Targets_Exp_Only"      )  

#######################
sink(file = paste(Sys.Date(),"_GSEA_AGE_session_Info.txt", sep =""))
sessionInfo()
sink()

