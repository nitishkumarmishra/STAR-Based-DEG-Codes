#!/bin/bash

source bashrc_rdna_variant_calling


bsub -q priority -n 10 -P ensmb -R "rusage[mem=25000]" -J star_ensmbl -oo STAR.ensmbl.out -eo STAR.ensmbl.err "sh generate_STAR_mRNA_genome_index_for_annotated_transcriptome.sh -p 10  mm39"

bsub -q priority -n 10 -P gencode -R "rusage[mem=25000]" -J star_gencode -oo STAR.gencode.out -eo STAR.gencode.err "sh generate_STAR_mRNA_genome_index_for_annotated_transcriptome.sh -p 10 GENCODEmm39"

