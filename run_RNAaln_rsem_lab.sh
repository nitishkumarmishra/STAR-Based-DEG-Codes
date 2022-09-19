#!/bin/bash

#src/rsem_eukaryotic_RNAseq_reads_to_annotated_genome.sh
## Modified by Nitish
#samples_n=0
for fq1 in fastq/*R1.fastq.gz; do

	#samples_n=$(( samples_n +1))
	#if test $samples_n -eq 1; then
	#	continue 
	#fi
        #echo -e "sample number :: $samples_n "

	fn=$(basename $fq1)
	dir=${fq1%%$fn}
	sample=${fn%%_R1.fastq.gz}
	fq2=$( ls ${dir}${sample}_R2.fastq.gz )
	echo $sample
	echo -e "[R1] $fq1"
	echo -e "[R2] $fq2"
	#PROG="$RPROF/src/rsem_eukaryotic_RNAseq_reads_to_annotated_genome.sh"
	PROG="/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/blancgrp_257956_RNAseq_total_stranded-1/rsem.test.sh"
	mkdir -p $sample/STAR
	#cmd="$PROG -1 $fq1 -2 $fq2 -l paired -strat RNAseq -mem 8 -o ${sample}/${sample}.STAR -t ${sample}.temp -p 4 -rsemTranBam 1"
	cmd="$PROG -1 $fq1 -2 $fq2 -l paired -mem 8 -rsemTranBam -starBam -stranded reverse -o ${sample}/STAR -t ${sample}.temp -p 4 -s mm39"
	echo -e "[CMD] $cmd"
	bsub -n 12 -J rna -e aln.err -o aln.out -q priority -P RSEM -R 'rusage[mem=25000]' $cmd	

done

#parameters:
#-u file	input single end fastq or fastq file of mixed mates for paired-end (will separate them appropriately). should be .gz
#-1 file	input fastq file for mate 1.  should be .gz
#-2 file	input fastq file for mate 2.  should be .gz
#-l layout	layout [paired, single].
#-s species	the species of the whole genome.  will load reference genome fasta and index if not provided (see "-ix" option).
#-strat  strategy	seuqencing strategy, e.g. RNAseq or Riboseq.
#-mem m	amount of memory on node, in GB.  JAVA with flag -XmxAG uses A+3 G memory (there is a secret 3 Gb overhead).  Will tailor accordingly to node memory. [default: assumes m = 8]

#If the following are provided, the output filename will be constructed and other parameters inferred.
#-runid r	runid


#-rsemTranBam 0/1	keep the transcriptome bam output from rsem [1] or discard [0]. Default: if strategy is ribosome profiling, then will keep, otherwise, discard.
#-take	will move the input fastq files to internal temporary storage.
#-replace	will replace fastq files to internal temporary storage, (if -take flag provided).
#-o file	output file prefix.  writes to <out>.bam
#-t pref	temporary file prefix.
#-p	number of threads for multithreading [default: 1].
#-k	keep intermediate files.



