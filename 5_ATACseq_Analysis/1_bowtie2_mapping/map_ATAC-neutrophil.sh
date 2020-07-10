WorkDir="/mnt/data/Data/FASTQ/ATACseq/Neutrophil_aging"

cd $WorkDir;

for f in $(find . -name '*_NGmerge_1.fastq.gz')
do
    f2=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/_NGmerge_2\.fastq\.gz/g');
    inf2="${WorkDir}/${f2}"

    of=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/\.bam/g');
    oFname="${WorkDir}/${of}"

    oflog=$(basename "${f}" | sed 's/_NGmerge_1\.fastq\.gz/\.map_log\.txt/g');
    LogFilename="${WorkDir}/${oflog}"

    bowtie2 --very-sensitive -k 10 -p 10 -X 2500 -x mm10 -1 $f -2 $inf2 --met-file $LogFilename | samtools view -b -S - > $oFname

done

