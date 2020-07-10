setwd('/Volumes/BB_Home/Neutrophils/Neutrophil_analysis/RNA-seq/GSEA')
options(stringsAsFactors = F)

# prep gmts for fast loading for GSEA with phenotest
# add sumamrized TF data
# filter TF data on expression in Neutrophils according to RNA-seq

library(DESeq2)
library(phenoTest)
library(qusage)

######### Process GMT files for use with ClsuterProfiler
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


convert_to_mouse <- function (my.geneset.list) {
  for (i in 1:length(my.geneset.list)) {
    my.geneset.list[[i]] <- firstup(tolower( my.geneset.list[[i]]  )) 
  }
  
  return(my.geneset.list)
}

#########
Sym.c2.cp.biocarta     <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c2.cp.biocarta.v7.0.symbols.gmt"                                              ) )
Sym.c2.cp.kegg         <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c2.cp.kegg.v7.0.symbols.gmt"                                                  ) )
Sym.c2.cp.reactome     <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c2.cp.reactome.v7.0.symbols.gmt"                                              ) )
Sym.c2.cp.pid          <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c2.cp.pid.v7.0.symbols.gmt"                                                   ) )
Sym.c2.cp.cgp          <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c2.cgp.v7.0.symbols.gmt"                                                      ) )
Sym.c2.cp              <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c2.cp.v7.0.symbols.gmt"                                                       ) )
Sym.c3.tft             <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c3.tft.v7.0.symbols.gmt"                                                      ) )
Sym.c7.immune          <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/c7.all.v7.0.symbols.gmt"                                                      ) )
Sym.h.all              <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/h.all.v7.0.symbols.gmt"                                                       ) )
Sym.hm.Chea            <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/CHEA_Harmonizome_2020-03-25.gmt"                                              ) )
Sym.hm.MPO             <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/Mammalian_Phenotype_Ontology_Harmonizome_2020-03-25.gmt"                      ) )
Sym.hm.HMDB            <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/HMDB_Harmonizome_2020-03-25.gmt"                                              ) )
Sym.hm.Interpro        <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/Interpro_Harmonizome_2020-03-25.gmt"                                          ) )
Sym.hm.LOCATE          <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/LOCATE_Harmonizome_2020-03-25.gmt"                                            ) )
Sym.hm.Panther         <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/PANTHER_Harmonizome_2020-03-25.gmt"                                           ) )
Sym.hm.PathComms       <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/Pathway_Commons_Interactions_Harmonizome_2020-03-25.gmt"                      ) )
Sym.hm.Wiki            <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/Wikipathways_Harmonizome_2020-03-25.gmt"                                      ) )
Sym.hm.dbGap           <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/dbGAP_Harmonizome_2020-03-25.gmt"                                             ) )
Sym.hm.GEO.tf.dwn      <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/GEO_TF_Perturb_DWN_Harmonizome_2020-03-25.gmt"                                ) )
Sym.hm.GEO.tf.up       <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/GEO_TF_Perturb_UP_Harmonizome_2020-03-25.gmt"                                 ) )
Sym.hm.Compartments    <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/COMPARTMENTS_Curated_Protein_Localization_Evidence_Harmonizome_2020-03-25.gmt") )
Sym.hm.kegg            <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/KEGG_Harmonizome_2020-03-25.gmt"                                              ) )
Sym.hm.EncodeTF        <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/ENCODE_TFs_Harmonizome_2020-03-25.gmt"                                        ) )
Sym.hm.JASPAR          <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/JASPAR_TF_targets_Harmonizome_2020-03-31.gmt"                                 ) )
Sym.GEO                <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/Pathways/TF-LOF_Expression_from_GEO_WITH_FOXO_2016-08-02_TF_targets.gmt"               ) )
Sym.miRTarBase         <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/miRNA/2020-04-09_mouse_miRTarBase7.0.gmt"                                              ) )
Sym.TargetScan         <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/miRNA/2020-04-10_mouse_TargetScan7_1.gmt"                                              ) )
Sym.ENS.GO.ALL         <- read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/ENSEMBL/2020-04-10_mouse_Ens99_GO_ALL.gmt"                                                               ) 
Sym.ENS.GO.BP          <- read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/ENSEMBL/2020-04-10_mouse_Ens99_GO_BP.gmt"                                                                ) 
Sym.ENS.GO.MF          <- read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/ENSEMBL/2020-04-10_mouse_Ens99_GO_MF.gmt"                                                                ) 
Sym.ENS.GO.CC          <- read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/ENSEMBL/2020-04-10_mouse_Ens99_GO_CC.gmt"                                                                ) 

Sym.hm.wiki.mouse      <- Sym.hm.Wiki[grep("musculus", names(Sym.hm.Wiki))]

Sym.TF.target.summary  <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/TF_GMT_cleanup/2020-05-20_TF_targets_ChEA_ENCODE_GEO_JASPAR_summary.gmt"              ))

save(Sym.c2.cp.biocarta   ,
     Sym.c2.cp.kegg       ,
     Sym.c2.cp.reactome   ,
     Sym.c2.cp.pid        ,
     Sym.c2.cp.cgp        ,
     Sym.c2.cp            ,
     Sym.c3.tft           ,
     Sym.c7.immune        ,
     Sym.h.all            ,
     Sym.hm.Chea          ,
     Sym.hm.MPO           ,
     Sym.hm.HMDB          ,
     Sym.hm.Interpro      ,
     Sym.hm.LOCATE        ,
     Sym.hm.Panther       ,
     Sym.hm.PathComms     ,
     Sym.hm.wiki.mouse    ,
     Sym.hm.dbGap         ,
     Sym.hm.GEO.tf.dwn    ,
     Sym.hm.GEO.tf.up     ,
     Sym.hm.Compartments  ,
     Sym.hm.kegg          ,
     Sym.hm.EncodeTF      ,
     Sym.hm.JASPAR        ,
     Sym.GEO              ,
     Sym.miRTarBase       ,
     Sym.TargetScan       ,
     Sym.ENS.GO.ALL       ,
     Sym.ENS.GO.BP        ,
     Sym.ENS.GO.MF        ,
     Sym.ENS.GO.CC        ,
     Sym.TF.target.summary,
     file = paste0(Sys.Date(),"_GeneSetCollections_for_Phenotest_GSEA.RData"))


############################################################################################
# Filter TF on expression

# DEseq2 counts filtered on expression
my.neutrophil.RNA.seq <- read.table("../DEseq2/STAR/2020-05-21_Neutrophils_DESeq2_Analysis_STAR_post_SVA__log2_counts_matrix_DEseq2_SVA.txt", sep = "\t", header = T)
my.neutrophil.exp     <- rownames(my.neutrophil.RNA.seq)

# read the TF summary
Sym.TF.target.summary  <- convert_to_mouse( read.gmt("/Volumes/BB_Home/PATHWAY_ANNOT/TF_GMT_cleanup/2020-05-20_TF_targets_ChEA_ENCODE_GEO_JASPAR_summary.gmt"              ))

# get names formatted for mouse
names(Sym.TF.target.summary) <- firstup(tolower(names(Sym.TF.target.summary))) # 391
 
my.exp.tfs <- names(Sym.TF.target.summary) %in% my.neutrophil.exp   # 293

Sym.TF.target.summary.exp <- Sym.TF.target.summary[my.exp.tfs]

save(Sym.TF.target.summary.exp   ,
     file = paste0(Sys.Date(),"Expressed_TF_summary_GeneSetCollections_for_Phenotest_GSEA.RData"))

