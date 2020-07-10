BED2GFF.pl -bed /Volumes/BB_Home/Neutrophils/Neutrophil_analysis/ATAC-seq_analysis/Peak_Calls/HOMER_peaks/Neutrophil_ATACseq_Peaks_meta.QCnoMt.bed -genome mm10
GFF_ATAC.FragSizes.pl -gff Neutrophil_ATACseq_Peaks_meta.QCnoMt.gff -file bam_for_frag_counting_v2.txt
GFF2TXT.pl -gff Neutrophil_ATACseq_Peaks_meta.QCnoMt.FragSizes.ANNOTATED.gff
