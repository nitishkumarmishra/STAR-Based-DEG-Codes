#!/bin/bash
#BSUB -P Bowtie2
#BSUB -oo Bowtie2.out -eo Bowtie2.err
#BSUB -n 5
#BSUB -M 10G
#BSUB -J bowtie3


source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv

bowtie2-build Sequence_rRNA_snRNA_tRNA_Duplicate.fa noncoding_RNA
mkdir -p Bowtie2/noncoding_RNA
mv noncoding_RNA* Bowtie2/noncoding_RNA/


