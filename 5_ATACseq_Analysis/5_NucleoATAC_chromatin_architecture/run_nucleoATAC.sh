nucleoatac run --cores 12 --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --bam Neutrophil_ATACseq_4m_M_DS_Merged.nucleoatac.bam  --fasta mm10.fa --write_all --out Neutrophil_ATACseq_4m_M_NucleoATAC
nucleoatac run --cores 12 --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --bam Neutrophil_ATACseq_4m_F_DS_Merged.nucleoatac.bam  --fasta mm10.fa --write_all --out Neutrophil_ATACseq_4m_F_NucleoATAC
nucleoatac run --cores 12 --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --bam Neutrophil_ATACseq_21m_M_DS_Merged.nucleoatac.bam --fasta mm10.fa --write_all --out Neutrophil_ATACseq_21m_M_NucleoATAC
nucleoatac run --cores 12 --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --bam Neutrophil_ATACseq_21m_F_DS_Merged.nucleoatac.bam --fasta mm10.fa --write_all --out Neutrophil_ATACseq_21m_F_NucleoATAC

nucleoatac merge --occpeaks Neutrophil_ATACseq_4m_M_NucleoATAC.occpeaks.bed.gz  --nucpos Neutrophil_ATACseq_4m_M_NucleoATAC.nucpos.bed.gz
nucleoatac merge --occpeaks Neutrophil_ATACseq_4m_F_NucleoATAC.occpeaks.bed.gz  --nucpos Neutrophil_ATACseq_4m_F_NucleoATAC.nucpos.bed.gz
nucleoatac merge --occpeaks Neutrophil_ATACseq_21m_M_NucleoATAC.occpeaks.bed.gz --nucpos Neutrophil_ATACseq_21m_M_NucleoATAC.nucpos.bed.gz
nucleoatac merge --occpeaks Neutrophil_ATACseq_21m_F_NucleoATAC.occpeaks.bed.gz --nucpos Neutrophil_ATACseq_21m_F_NucleoATAC.nucpos.bed.gz

nucleoatac nfr --bam Neutrophil_ATACseq_4m_M_DS_Merged.nucleoatac.bam  --fasta mm10.fa --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --occ_track Neutrophil_ATACseq_4m_M_NucleoATAC.occ.bedgraph.gz  --calls Neutrophil_ATACseq_4m_M_NucleoATAC.nucmap_combined.bed.gz --cores 12 --out Neutrophil_ATACseq_4m_M_NucleoATAC
nucleoatac nfr --bam Neutrophil_ATACseq_4m_F_DS_Merged.nucleoatac.bam  --fasta mm10.fa --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --occ_track Neutrophil_ATACseq_4m_F_NucleoATAC.occ.bedgraph.gz  --calls Neutrophil_ATACseq_4m_F_NucleoATAC.nucmap_combined.bed.gz --cores 12 --out Neutrophil_ATACseq_4m_F_NucleoATAC
nucleoatac nfr --bam Neutrophil_ATACseq_21m_M_DS_Merged.nucleoatac.bam --fasta mm10.fa --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --occ_track Neutrophil_ATACseq_21m_M_NucleoATAC.occ.bedgraph.gz --calls Neutrophil_ATACseq_21m_M_NucleoATAC.nucmap_combined.bed.gz --cores 12 --out Neutrophil_ATACseq_21m_M_NucleoATAC
nucleoatac nfr --bam Neutrophil_ATACseq_21m_F_DS_Merged.nucleoatac.bam --fasta mm10.fa --bed Neutrophil_ATACseq_regions_for_NucleoATAC.bed --occ_track Neutrophil_ATACseq_21m_F_NucleoATAC.occ.bedgraph.gz --calls Neutrophil_ATACseq_21m_F_NucleoATAC.nucmap_combined.bed.gz --cores 12 --out Neutrophil_ATACseq_21m_F_NucleoATAC
 