#!/bin/bash

if [[ "$#" -lt 3 ]]
then
    echo "$(basename $0) [BamDir] [oDir] "  1>&2
    echo "   [Ass]: assembly" 1>&2 
    echo "   [BamDir]: folder with bam files" 1>&2
    echo "   [oDir]: output directory" 1>&2

    exit 1
fi

Ass=$(echo $1 | sed 's:/$::g')
BamDir=$(echo $2 | sed 's:/$::g')
oDir=$(echo $3 | sed 's:/$::g')

# make output directory if it doesnt exist
[[ ! -d "${oDir}" ]] && mkdir "${oDir}"
[[ ! -d "${oDir}/BAM_clean" ]] && mkdir "${oDir}/BAM_clean"

##########################################
# 1. convert SAM into BAM, clean up alignemnts
for f in $(find "${BamDir}" -name '*.bam')
do
	filePath="${oDir}/BAM_clean"
	
	# filter to retain mapped reads
	fileName=$(basename "${f}" | sed 's/\.bam/\.mp\.bam/g');
    oFname="${filePath}/${fileName}"
    samtools view -b -F 4 -q 15 $f > $oFname # mapped only
    
    # sort mapped reads (space gain)
	outBamPre=$(basename "${f}" | sed 's/\.bam/\.mp\.srt/g');
	oFname2b="${filePath}/${outBamPre}"
	outBam2=$(basename "${f}" | sed 's/\.bam/\.mp\.srt\.bam/g');
	oFname2="${filePath}/${outBam2}"

	#echo $oFname2b
	#echo $oFname2
	
	samtools sort -T $oFname2b -o $oFname2 $oFname 
		
	# clean up pre-sort file
	rm $oFname
	
	# index sorted bam file
	samtools index $oFname2

	# These flags exclude unmapped reads, secondary alignments,
	# and QC failures.
	outFile3=$(basename "${f}" | sed 's/\.bam/\.rmdup\.bam/g');
	oFname3="${filePath}/${outFile3}"
	
	samtools rmdup $oFname2 $oFname3
	
	# determine duplication rates and store in report file
	getDupRate_PE.pl $oFname2 $oFname3 >> $filePath/Libraries_PCR_duplication_rates_report.txt

	
	outHOM=$(basename "${f}" | sed 's/\.bam/_Filtered_TAGs/g');
	oDir4="${filePath}/${outHOM}"
	
	makeTagDirectory $oDir4 $oFname3 -genome $Ass -format sam -keepOne
	
    
done
echo "Finished bam clean up"

