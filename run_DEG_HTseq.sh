#!/bin/bash
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv
module load python/3.7.0
module load R/4.1.0
module load gcc/9.1.0


bsub -n 5 -J rna -e DEG_HTseq.err -o DEG_HTseq.out -q priority -P RSEM -R 'rusage[mem=5000]' "Rscript run_geneDE_HTSeq.R"
