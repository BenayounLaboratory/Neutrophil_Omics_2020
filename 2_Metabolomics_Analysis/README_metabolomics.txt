README -- metabolomics analysis
############################

	1. Data processing from all modes and columns (PQI_data_processing_RPLC_pos.R, PQI_data_processing_RPLC_neg.R, PQI_data_processing_HILIC_pos.R, PQI_data_processing_HILIC_neg.R, Merge datasets.R)
	2. Analyze data with limma, normalize with BCA readings (Process_metabolomics_limma_v2.R, Functions_metabolomics_limma_v2.R)
	3. Use MetaboanalystR with mummichog algorithm
			- get PSEA information (run_metaboanalyst_R_v2.R)
			- Make bubble chart ouput (Get_MetabBubble_chart_v2.R)
	4. Integrated pathway analysis with IMPaLA
			- online analysis
			- parsing in R (parse_plot_Impala.R)
