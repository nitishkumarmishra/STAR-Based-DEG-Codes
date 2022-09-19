#!/bin/bash
#BSUB -P Bowtie1_pep
#BSUB -oo Bowtie1_pep.out -eo Bowtie1_pep.err
#BSUB -n 5
#BSUB -M 10G
#BSUB -J tRNA

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



bowtie-build Sequence_rRNA_snRNA_tRNA_Duplicate.fa noncoding_RNA
mkdir -p Bowtie1.3.0/noncoding_RNA
mv noncoding_RNA* Bowtie1.3.0/noncoding_RNA/


## -f 1 is transcrpt and -f 2 for geneID. If I am using -f 2 then sed 's/ENS/>ENS/, because 1st part have > but not second.
## First is transcript ID and second is geneID
cut -d '|' -f 1 gencode.vM27.pc_transcripts.fa |sed 's/ENS/ENS/g' > gencode.vM27.pc_transcripts_new.fa
#cut -d '|' -f 2 gencode.vM27.pc_transcripts.fa |sed 's/ENS/>ENS/g' > gencode.vM27.pc_transcripts_new.fa
bowtie-build --offrate gencode.vM27.pc_transcripts_new.fa gencode.vM27.pc
mkdir -p Bowtie1.3.0/gencode.vM27.pc
mv gencode.vM27.pc* Bowtie1.3.0/gencode.vM27.pc/
#bowtie -q -p 18 -v 2 mm39-mature-tRNA infile.fastq --un unaligned_output.fastq --al aligned_output.fastq
