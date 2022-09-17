#!/bin/sh
#SBATCH --time=95:00:00          # Run time in hh:mm:ss
#SBATCH --mem=50g       # Maximum memory required per CPU (in megabytes)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --partition=batch
#SBATCH --job-name=HTSeq_1
#SBATCH --output=slurmout/Job_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurmout/Job_%A_%a.err # File to which STDERR will be written

module load compiler/gcc
module load star/2.7
module load htseq/py37/0.11

sh HTSeq_1.sh
