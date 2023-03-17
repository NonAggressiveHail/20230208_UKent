#!/bin/bash

#SBATCH --job-name=p_aeruginosa_pangenomics
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=14000
#SBATCH --error=p_aeruginosa_pangenomics.err
#SBATCH --out=p_aeruginosa_pangenomics.out
#SBATCH --time=84:00:00

# conda needs special stuff to start via slurm
eval "$(/shared/home/jgh8/miniconda3/bin/conda shell.bash hook)"

conda activate anvio-7.1

# Main program
bash 20230309_anvio_pangenomics.sh -g ../data/RefSeq_pseudomonas -n $SLURM_CPUS_PER_TASK -o ../data/test_anvio
