samtools view Neutrophil_ATACseq_4m_M_DS_Merged.srt.bam  | grep -v 'chrM' | samtools view -SbT /Users/berenice/Softwares/Genomes/mm10.fa -o Neutrophil_ATACseq_4m_M_DS_Merged.srtnoMT.bam  -
samtools view Neutrophil_ATACseq_4m_F_DS_Merged.srt.bam  | grep -v 'chrM' | samtools view -SbT /Users/berenice/Softwares/Genomes/mm10.fa -o Neutrophil_ATACseq_4m_F_DS_Merged.srtnoMT.bam  -
samtools view Neutrophil_ATACseq_21m_M_DS_Merged.srt.bam | grep -v 'chrM' | samtools view -SbT /Users/berenice/Softwares/Genomes/mm10.fa -o Neutrophil_ATACseq_21m_M_DS_Merged.srtnoMT.bam -
samtools view Neutrophil_ATACseq_21m_F_DS_Merged.srt.bam | grep -v 'chrM' | samtools view -SbT /Users/berenice/Softwares/Genomes/mm10.fa -o Neutrophil_ATACseq_21m_F_DS_Merged.srtnoMT.bam -

# get BED with the recommended 250bp slop
slopBed -i /Volumes/BB_Home/Neutrophils/Neutrophil_analysis/ATAC-seq_analysis/Peak_Calls/HOMER_peaks/Neutrophil_ATACseq_Peaks_meta.QCnoMt.bed -g /Users/berenice/Softwares/bedtools-master/genomes/mouse.mm10.genome -b 250 > Neutrophil_ATACseq_Peaks_meta.QCnoMt.SLOP.bed

# also extract TSS coordinates from expressed genes (-750,+750)

annotatePeaks.pl tss mm10 -list 2020-04-06_Neutrophil_Expressed_Genes_from_STAR.txt -size -750,750 -noann > HOMER_Neutrophil_Expressed_TSS.txt

cat HOMER_Neutrophil_Expressed_TSS.txt | tail -12003 | cut -f 2,3,4,16 | sortBed -i - > HOMER_Neutrophil_Expressed_TSS_coords.bed

mergeBed -i HOMER_Neutrophil_Expressed_TSS_coords.bed > HOMER_Neutrophil_Expressed_TSS_coords.QC.bed

cat Neutrophil_ATACseq_Peaks_meta.QCnoMt.SLOP.bed HOMER_Neutrophil_Expressed_TSS_coords.QC.bed | wc -l

cat Neutrophil_ATACseq_Peaks_meta.QCnoMt.SLOP.bed HOMER_Neutrophil_Expressed_TSS_coords.QC.bed | sortBed -i - | mergeBed -d 50 -i - > Neutrophil_ATACseq_regions_for_NucleoATAC.bed