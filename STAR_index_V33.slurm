#!/bin/sh
#SBATCH --ntasks=16 # Number of task
#SBATCH --time=30:00:00          # Run time in hh:mm:ss
#SBATCH --mem=50g       # Maximum memory required per CPU (in megabytes)
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --job-name=STAR_hg38
#SBATCH --error=/work/unmc_gudalab/nitish123/STAR_hg38/job.%J.err
#SBATCH --output=/work/unmc_gudalab/nitish123/STAR_hg38/job.%J.out

module load star/2.7
STAR --runMode genomeGenerate --genomeDir /work/unmc_gudalab/nitish123/STAR_hg38/STAR_hg38_V33 --genomeFastaFiles /work/unmc_gudalab/nitish123/STAR_hg38/hg38.fa --sjdbGTFfile /work/unmc_gudalab/nitish123/STAR_hg38/gencode.v33.annotation.gtf --sjdbOverhang 100 --runThreadN 8
