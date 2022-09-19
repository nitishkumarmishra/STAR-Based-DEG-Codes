#!/bin/bash
#BSUB -P Trimgalore
#BSUB -oo Trimgalore.out -eo Trimgalore.err
#BSUB -n 15
#BSUB -M 20G
#BSUB -J Trimgalore

NUMBER_THREAD=10
####################################################
############### Parsing FASTQ file #################
PREFIX=$1
FQ1=`ls|grep $PREFIX | grep "_R1_"| grep "fastq.gz"`
FQ2=`ls|grep $PREFIX | grep "_R2_"| grep "fastq.gz"`

####################################################
######### Trimgalore adapter trimming ##############
### Conda environ for latest cutadapt and trmgalore.
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv
mkdir -p ./FASTQC/TRIM #Output file for FASTQ
mkdir -p ./TRIM ##
OUTPUT_DIR=./TRIM

#/home/nmishra/TrimGalore-0.6.6/trim_galore \
trim_galore \
-j $NUMBER_THREAD \
--clip_R1 10 \
--clip_R2 10 \
--basename $PREFIX \
--output_dir ${OUTPUT_DIR} \
--paired \
--fastqc_args "--outdir ./FASTQC/TRIM" $FQ1 $FQ2

