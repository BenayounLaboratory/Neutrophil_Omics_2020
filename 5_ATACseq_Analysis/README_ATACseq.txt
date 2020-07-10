README -- ATAC-seq analysis
############################

	1. Reads trimmed with NGMerge and mapped to mm10 with bowtie2 
			- map with bowtie2 (map_ATAC-neutrophil.sh)
			- remove duplicates  (clean_BAM_PE.sh)
			- downsample so all libraries are "equal" (downsample_files_for_comparison.sh,Down_sampling_bam_file_rev.pl)
			- create HOMER TAG directories (make_homer_dirs.sh, get_merged_TAG_dirs.sh)
	2. Call peaks with HOMER (metapeaks over collapsed biological replicates, then common to each condition)
			- call_ATAC-peaks.sh
	3. Run differential accessibility analysis with Diffbind/DEseq2
			- get normalized data using diffbind (Diffbind_Analyze_ATAC_v2.R)
			- get differentially accessible data using DEseq2 (process_sex_ATACseq_datasets_v4.R)
	4. Plot genomic coordinates of sex-dimorphic genes (Omicscircos R package)
			- plot positions (Plos_circle_plot_ATAC.R)
	5. Chromatin architecture with NucleoATAC
			- prep BAM files and input bed according to developer instructions (prep_BAM_files.sh)
			- run nucleoATAC pipeline (run_nucleoATAC.sh)
			- median TSS occupancy signal:
				+ Extract signal from wig with HOMER: 
					annotatePeaks.pl tss mm10 -list 2020-04-06_Neutrophil_Expressed_Genes_from_STAR.txt -size -500,500 -hist 10 -ghist -bedGraph Neutrophil_ATACseq_4m_F_NucleoATAC.occ.bedgraph Neutrophil_ATACseq_4m_M_NucleoATAC.occ.bedgraph Neutrophil_ATACseq_21m_F_NucleoATAC.occ.bedgraph Neutrophil_ATACseq_21m_M_NucleoATAC.occ.bedgraph > HOMER_Neutrophil_NucleoATAC_occupancy_4mFM_21mFM.txt
				+ plot median of occupancy (Plot_median_occupancy.R)
			- boxplot occupancy, fuzziness, etc from NucleoATAC output (Process_NucleoATAC_results_v2.R)
			- Fragment lengths captured with gfftools (Analyze_nucleoATAC_output.sh, Analyze_fragments_sizes_gfftools.R)