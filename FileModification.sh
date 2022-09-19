#!/bin/bash
#BSUB -P Edit_file
#BSUB -oo Edit_file.out -eo Edit_file.err
#BSUB -n 5
#BSUB -M 10G
#BSUB -J Edit_file

#module load bowtie/1.2
#module load fastx_toolkit/0.0.13
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv


fasta_formatter -w 0 -i Mus_musculus.GRCm39.ncrna.fa -o mature_idx.fa
#sed -i -e  '/^>/s/ .*//' mature_idx.fa
fasta_formatter -t  -i mature_idx.fa -o mature_idx_1.fa ## Convert file in single line FASTA format, alternate awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'
sed -i -e 's/^/>/g' mature_idx_1.fa
grep -i -E "snRNA|rRNA" mature_idx_1.fa >mature_idx_2.fa
sed -i -e 's/\t/\n/g' mature_idx_2.fa
cut -d ' ' -f1 mature_idx_2.fa > mature_idx_3.fa
mv mature_idx_3.fa Mus_musculus.GRCm39_tRNA_snRNA.fa
#rm mature_idx.fa mature_idx_1.fa mature_idx_2.fa mature_idx_3.fa


seqkit seq --rna2dna mm39-mature-tRNAs.fa > mature_igenome.fa
fasta_formatter -w 0 -i mature_igenome.fa -o mature_idx.fa
sed -i -e  '/^>/s/ .*//' mature_idx.fa
fasta_formatter -t  -i mature_idx.fa -o mature_idx_1.fa
sed -i -e 's/CCA$//g' mature_idx_1.fa 
sed -i -e 's/$/CCA/g' mature_idx_1.fa
sed -i -e 's/^/>/g' mature_idx_1.fa
sed -i -e 's/\t/\n/g' mature_idx_1.fa
mv mature_idx_1.fa mm39-mature-tRNAs_3P_CCA.fa
cat Mus_musculus.GRCm39_tRNA_snRNA.fa mm39-mature-tRNAs_3P_CCA.fa >Sequence_rRNA_snRNA_tRNA.fa
seqkit rmdup -s Sequence_rRNA_snRNA_tRNA.fa > Sequence_rRNA_snRNA_tRNA_Duplicate.fa ## Remove duplicated sequence


rm mature_idx.fa mature_idx_1.fa mature_idx_2.fa mature_idx_3.fa


bowtie-build Sequence_rRNA_snRNA_tRNA_Duplicate.fa noncoding_RNA
mkdir -p Bowtie1.3.0/noncoding_RNA
mv noncoding_RNA* Bowtie1.3.0/noncoding_RNA/


INDEX=/research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/GeneSeq/tRNA/Bowtie2.4.4/mm39-tRNA/
### Conda environ for latest bowtie2 ########.
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv

##### Bowtie2 alignment against tRNA (not mature)
bowtie2 --quiet --min-score G,1,8 --local -D 20 -R 3 -N 1 -L 10 -i S,1,0.5 -p 10 -x ${INDEX}immature_tRNA  -1 2113859_JBQS011_val_1.fq.gz -2 2113859_JBQS011_val_2.fq.gz -S 2113859.sam

## -p==paired end; -M==Multimapping; O==AllowMuliOverlap; T==Number of thread; S==Strand Specific (0=none, 1=forward, 2=reverse)
featureCounts -p -M -O -Q 10 -S 2 -T 6 --fraction -a /research/groups/blancgrp/home/nmishra/EMT_data/Small_RNA/FASTQ/Trial/mm39-tRNAs_exclude_Pseudogene_High_Confidense.gtf -o 2113859.counts.txt 2113859.sam
~

