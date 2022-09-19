#!/bin/bash

#for PREFIX in $(ls |grep "fastq")
for PREFIX in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | cut -d _ -f 1,2,3,4  ; done | sort | uniq)
do
	#echo "${PREFIX}"
	#echo -e "\n"
	echo "****************** Sample $PREFIX submitted *********************"
        #echo $PREFIX
        bsub -q priority -n 14 -P Auto_trim -R "rusage[mem=25000]" -J Trim_fq -oo %J.out -eo %J.err "sh Trimgalore.sh $PREFIX"
        echo "*************** Gene counts for sample $PREFIX ******************"
        echo -e "#########################################################################""\n"

done

