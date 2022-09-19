#!/bin/bash
#BSUB -P STAR_COUNT
#BSUB -oo STAR_COUNT.out -eo STAR_COUNT.err
#BSUB -n 15
#BSUB -M 30G
#BSUB -J STAR_COUNT

####################################################
############ GTF and STAR genome Index #############
GENOME=/home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/STAR-index/2.7.9a
GTF=/home/nmishra/REF/Mus_musculus/REF_GRCm39_vM27/REF/gencode.vM27.annotation.gtf
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


TRIM1=${PREFIX}_val_1.fq.gz
TRIM2=${PREFIX}_val_2.fq.gz
####################################################
############## STAR GeneCounts #####################
STAR \
--genomeDir $GENOME \
--readFilesCommand zcat \
--runThreadN $NUMBER_THREAD \
--sjdbGTFfile $GTF \
--alignEndsType EndToEnd \
--alignSJoverhangMin 8 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--alignSJDBoverhangMin 10 \
--alignIntronMin 20 \
--alignIntronMax 100000 \
--alignMatesGapMax 100000 \
--outFileNamePrefix ${PREFIX}_ \
--outFilterMismatchNmax 50 \
--outSAMattributes NH HI AS nM NM ch \
--outFilterMultimapScoreRange 10 \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 999 \
--outFilterMismatchNoverReadLmax 0.33 \
--outSAMunmapped Within \
--outFilterMismatchNoverLmax 0.1 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--alignSJDBoverhangMin 1000 \
--outWigType wiggle \
--outWigStrand Stranded \
--outWigNorm RPM \
--quantMode TranscriptomeSAM GeneCounts \
--twopassMode Basic \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 12 \
--chimSegmentReadGapMax 3 \
--chimMultimapNmax 10 \
--chimMultimapScoreRange 10 \
--chimNonchimScoreDropMin 10 \
--chimOutJunctionFormat 1 \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1 \
--readFilesIn $TRIM1 $TRIM2

####################################################
############## RSEM GeneCounts #####################
# $1 == STAR alignments bam file
# $2 == name of STAR index
# $3 == sample name
mkdir -p ./RSEM_COUNTS/
RSEM_OUT=./RSEM_COUNTS

rsem-calculate-expression \
-p 8 \
--no-bam-output  \
--alignments  \
--paired-end  \
--strandedness reverse \
${PREFIX}_Aligned.toTranscriptome.out.bam \
/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/GENOME_data/rsem/tRNA/vM27_tRNA \
${PREFIX}.RSEM

