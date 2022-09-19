#!/bin/bash
#BSUB -P fastqc
#BSUB -oo fastqc.out -eo fastqc.err
#BSUB -n 5
#BSUB -M 15G
#BSUB -J fastqc

mkdir -p fastqc_results/
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv

for i in *.fastq.gz
do
        fastqc -o fastqc_results/ $i
done

multiqc fastqc_results/

