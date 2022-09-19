#!/bin/bash

source bashrc_rdna_variant_calling


bsub -q priority -n 10 -P RSEM -R "rusage[mem=25000]" -J rsem_idx -oo %J.out -eo %J.err "sh generate_rsem_transcriptome_index.sh mm39"

bsub -q priority -n 10 -P RSEM -R "rusage[mem=25000]" -J rsem_idx -oo %J.out -eo %J.err "sh generate_rsem_transcriptome_index.sh GENCODEmm39"

