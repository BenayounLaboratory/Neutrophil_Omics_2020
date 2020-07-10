README  -    Impala analysis
############################

	1. Get FDR <5% genes from DEseq2
	2. Get manually validated metabolites at FDR <5% (HMBD IDs only)
	3. Get significant lipids at FDR <5%
		- use lipid annotation file to transfer HMDB IDs to results based on LION lipid name (Lists/get_lipid_Lists_HMDB.R)
		- only lipids with existing HMDB IDs (Kevin_Lipid_annotation.txt) are used
	4. Concatenate lists for metabolites and lipids as one input list of HMDBIDs
	5. Run - Impala version 12, access 2020-06-30
		- genes: HNGC/Symbols
		- metabolites/lipids:HMDB IDs
