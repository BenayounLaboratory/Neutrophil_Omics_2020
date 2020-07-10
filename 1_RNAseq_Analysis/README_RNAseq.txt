README -- RNA-seq analysis
############################

	1. Reads trimmed with trim_galore and mapped to mm10 with STAR (map_with_star.sh)
	2. Count reads mapping to genes using feature counts (count_reads.sh) and prepare TAG directory for HOMER/repeat analysis (get_RNA_TAG_directories.sh)
	3. Run DESeq2 analysis for DE genes as a function of age and/or sex (process_sex_RNAseq_datasets_STAR_clean)
	4. Run GSEA enrichment analysis:
			- prep gmt gene sets for GSEA (Prep_and_Load_Gene_Sets.R)
			- run GSEA (Run_GSEA_analysis_Aging_v2.R,Run_GSEA_analysis_Sex_v2.R,Functions_GSEA_v2.R)
			- plot top 10 each direction using bubble chart (Get_Bubble_charts_v3.R, Plot_bubble_chart_functions_v4.R)
	5. Gene expression heatmaps (Plot_Heatmap.R)
	6. Plot genomic coordinates of sex-dimorphic genes (Omicscircos R package)
			- get mm10 TSS positions from HOMER
			- plot positions (Plot_circle_plot.R)
	7. Repeat expression analysis using HOMER
			- get counts in repeat and transcript (for normalization) from HOMER (Count_repeats.sh)
			- run DEseq2 pipeline to identify significant repeats (process_sex_Repeat_RNA.R)