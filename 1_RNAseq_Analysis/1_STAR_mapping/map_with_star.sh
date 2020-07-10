WorkDir="/mnt/data/Data/FASTQ"

cd $WorkDir;

export GENIN="/home/benayoun/Softwares/Genomes/MM10"

for f in $(find . -name '*R1_val_1.fq.gz')
do
    f2=$(basename "${f}" | sed 's/R1_val_1\.fq\.gz/R2_val_2\.fq\.gz/g');
    inf2="${WorkDir}/OB/${f2}"
    of=$(basename "${f}" | sed 's/R1_val_1\.fq\.gz//g');
    oFname="${WorkDir}/${of}"
	STAR --genomeDir $GENIN --readFilesIn $f $inf2 --readFilesCommand zcat --runThreadN 6 --outFilterMultimapNmax 200 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $oFname
done
