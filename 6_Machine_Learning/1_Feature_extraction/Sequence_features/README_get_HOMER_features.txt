annotatePeaks.pl tss mm10 -list 2020-04-06_Neutrophil_Expressed_Genes_from_STAR.txt -size -500,150 -CpG > 2020-04-28_HOMER_Neutrophil_Promoter_features.txt
annotatePeaks.pl rna mm10 -list 2020-04-06_Neutrophil_Expressed_Genes_from_STAR.txt -size given -CpG > 2020-04-28_HOMER_Neutrophil_Gene_features.txt
