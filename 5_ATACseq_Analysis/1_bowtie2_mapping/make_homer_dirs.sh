for f in $(find "." -name '*.bam')
do
    samtools index $f
done


makeTagDirectory Neutrophil_ATACseq_4m_F_101_TAGs/ Neutrophil_ATACseq_4m_F_101_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_F_102_TAGs/ Neutrophil_ATACseq_4m_F_102_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_F_103_TAGs/ Neutrophil_ATACseq_4m_F_103_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_F_104_TAGs/ Neutrophil_ATACseq_4m_F_104_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_F_105_TAGs/ Neutrophil_ATACseq_4m_F_105_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_M_111_TAGs/ Neutrophil_ATACseq_4m_M_111_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_M_112_TAGs/ Neutrophil_ATACseq_4m_M_112_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_M_113_TAGs/ Neutrophil_ATACseq_4m_M_113_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_M_114_TAGs/ Neutrophil_ATACseq_4m_M_114_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_4m_M_115_TAGs/ Neutrophil_ATACseq_4m_M_115_MergedLanes.rmdup_DS_45716676.bam   -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_F_106_TAGs/ Neutrophil_ATACseq_21m_F_106_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_F_107_TAGs/ Neutrophil_ATACseq_21m_F_107_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_F_108_TAGs/ Neutrophil_ATACseq_21m_F_108_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_F_109_TAGs/ Neutrophil_ATACseq_21m_F_109_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_F_110_TAGs/ Neutrophil_ATACseq_21m_F_110_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_M_116_TAGs/ Neutrophil_ATACseq_21m_M_116_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_M_117_TAGs/ Neutrophil_ATACseq_21m_M_117_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_M_118_TAGs/ Neutrophil_ATACseq_21m_M_118_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_M_119_TAGs/ Neutrophil_ATACseq_21m_M_119_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne
makeTagDirectory Neutrophil_ATACseq_21m_M_120_TAGs/ Neutrophil_ATACseq_21m_M_120_MergedLanes.rmdup_DS_45716676.bam -genome mm10 -format sam -keepOne

makeTagDirectory Neutrophil_ATACseq_4m_F_ALL_TAGs/  -d Neutrophil_ATACseq_4m_F_101_TAGs/ Neutrophil_ATACseq_4m_F_102_TAGs/ Neutrophil_ATACseq_4m_F_103_TAGs/ Neutrophil_ATACseq_4m_F_104_TAGs/ Neutrophil_ATACseq_4m_F_105_TAGs/  
makeTagDirectory Neutrophil_ATACseq_20m_F_ALL_TAGs/ -d Neutrophil_ATACseq_21m_F_106_TAGs/ Neutrophil_ATACseq_21m_F_107_TAGs/ Neutrophil_ATACseq_21m_F_108_TAGs/ Neutrophil_ATACseq_21m_F_109_TAGs/ Neutrophil_ATACseq_21m_F_110_TAGs/ 
makeTagDirectory Neutrophil_ATACseq_4m_M_ALL_TAGs/  -d Neutrophil_ATACseq_4m_M_111_TAGs/ Neutrophil_ATACseq_4m_M_112_TAGs/ Neutrophil_ATACseq_4m_M_113_TAGs/ Neutrophil_ATACseq_4m_M_114_TAGs/ Neutrophil_ATACseq_4m_M_115_TAGs/  
makeTagDirectory Neutrophil_ATACseq_20m_M_ALL_TAGs/ -d Neutrophil_ATACseq_21m_M_116_TAGs/ Neutrophil_ATACseq_21m_M_117_TAGs/ Neutrophil_ATACseq_21m_M_118_TAGs/ Neutrophil_ATACseq_21m_M_119_TAGs/ Neutrophil_ATACseq_21m_M_120_TAGs/ 

makeUCSCfile Neutrophil_ATACseq_4m_F_ALL_TAGs/  -fragLength given -res 3 -avg -raw -o Neutrophil_ATACseq_4m_F_ALL.bedgraph  -style dnase -skipChr chrM
makeUCSCfile Neutrophil_ATACseq_20m_F_ALL_TAGs/ -fragLength given -res 3 -avg -raw -o Neutrophil_ATACseq_20m_F_ALL.bedgraph -style dnase -skipChr chrM
makeUCSCfile Neutrophil_ATACseq_4m_M_ALL_TAGs/  -fragLength given -res 3 -avg -raw -o Neutrophil_ATACseq_4m_M_ALL.bedgraph  -style dnase -skipChr chrM
makeUCSCfile Neutrophil_ATACseq_20m_M_ALL_TAGs/ -fragLength given -res 3 -avg -raw -o Neutrophil_ATACseq_20m_M_ALL.bedgraph -style dnase -skipChr chrM
  