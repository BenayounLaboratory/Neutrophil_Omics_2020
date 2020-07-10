README


######### TF-LOF
# get FOXOs from Genome Research gmt
#     - save into FOXO-LOF_GEO.gmt

# Concatenate GEO TF-LOF for processing
cat GEO_TF_Perturb_DWN_Harmonizome_2020-03-25.gmt GEO_TF_Perturb_UP_Harmonizome_2020-03-25.gmt FOXO-LOF_GEO.gmt > GEO_TF_Perturb_Harmonizome_2020-03-25_UP_DOWN_FOXO.gmt

# make summary GMT
./clean_GEO_gmt.pl GEO_TF_Perturb_Harmonizome_2020-03-25_UP_DOWN_FOXO.gmt