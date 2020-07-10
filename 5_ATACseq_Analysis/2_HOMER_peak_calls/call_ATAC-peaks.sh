findPeaks ../BOWTIE2_mapped/Neutrophil_ATACseq_4m_M_Merged_TAGs -style dnase  -o Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.pos  -region 
findPeaks ../BOWTIE2_mapped/Neutrophil_ATACseq_4m_F_Merged_TAGs -style dnase  -o Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.pos  -region 
findPeaks ../BOWTIE2_mapped/Neutrophil_ATACseq_21m_F_Merged_TAGs -style dnase -o Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.pos -region 
findPeaks ../BOWTIE2_mapped/Neutrophil_ATACseq_21m_M_Merged_TAGs -style dnase -o Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.pos -region 


pos2bed.pl -o Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.bed   Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.pos  
pos2bed.pl -o Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.bed   Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.pos  
pos2bed.pl -o Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.bed   Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.pos  
pos2bed.pl -o Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.bed   Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.pos  

#	Output File: Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.bed
#	Converted 81683 peaks total

#	Output File: Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.bed
#	Converted 85848 peaks total\

#	Output File: Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.bed
#	Converted 69861 peaks total

#	Output File: Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.bed
#	Converted 72517 peaks total

sortBed -i Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.bed  > Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.srt.bed  
sortBed -i Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.bed  > Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.srt.bed  
sortBed -i Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.bed > Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.srt.bed 
sortBed -i Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.bed > Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.srt.bed 

intersectBed -v -a Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.srt.bed  -b /Users/berenice/Softwares/Genomes/mm10-blacklist.v2.bed > Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.QC.bed 
intersectBed -v -a Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.srt.bed  -b /Users/berenice/Softwares/Genomes/mm10-blacklist.v2.bed > Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.QC.bed 
intersectBed -v -a Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.srt.bed -b /Users/berenice/Softwares/Genomes/mm10-blacklist.v2.bed > Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.QC.bed
intersectBed -v -a Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.srt.bed -b /Users/berenice/Softwares/Genomes/mm10-blacklist.v2.bed > Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.QC.bed

#multiIntersectBed -i Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.QC.bed Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.QC.bed Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.QC.bed Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.QC.bed \
#   | perl -lane 'print if $F[3]>=2' | cut -f 1,2,3 | mergeBed -d 50 -i - | perl -lane 'print "$_\tATAC_Peak_$."' > Neutrophil_ATACseq_Peaks_across_2_or_more_conditions.bed
#
#cat Neutrophil_ATACseq_Peaks_across_2_or_more_conditions.bed | cut -f 1,2,3 > Neutrophil_ATACseq_Peaks_across_2_or_more_conditions.NoName.bed
#

cat Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.QC.bed  | grep -v 'chrM' > Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.QCnoMt.bed 
cat Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.QC.bed  | grep -v 'chrM' > Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.QCnoMt.bed 
cat Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.QC.bed | grep -v 'chrM' > Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.QCnoMt.bed
cat Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.QC.bed | grep -v 'chrM' > Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.QCnoMt.bed

#68040 Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.QCnoMt.bed
#70640 Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.QCnoMt.bed
#83727 Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.QCnoMt.bed
#79742 Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.QCnoMt.bed

cat Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.QCnoMt.bed Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.QCnoMt.bed  Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.QCnoMt.bed Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.QCnoMt.bed | sortBed -i - | mergeBed -d 50 -i - > Neutrophil_ATACseq_Peaks_meta.QCnoMt.bed
# 95704

multiIntersectBed -i Neutrophil_ATACseq_4m_M_MERGED_homer_peaks.QCnoMt.bed Neutrophil_ATACseq_4m_F_MERGED_homer_peaks.QCnoMt.bed Neutrophil_ATACseq_21m_F_MERGED_homer_peaks.QCnoMt.bed Neutrophil_ATACseq_21m_M_MERGED_homer_peaks.QCnoMt.bed \
   | perl -lane 'print if $F[3]>=2' | cut -f 1,2,3 | mergeBed -d 50 -i - | perl -lane 'print "$_\tATAC_Peak_$."' > Neutrophil_ATACseq_Peaks_across_2_or_more_conditions.QC.bed
#   78148 Neutrophil_ATACseq_Peaks_across_2_or_more_conditions.QC.bed
