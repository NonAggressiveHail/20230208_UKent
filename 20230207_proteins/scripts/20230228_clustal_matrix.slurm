#!/bin/bash

#SBATCH --job-name=clustal
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --error=clustal.err
#SBATCH --out=clustal.out
#SBATCH --time=24:00:00

# conda needs special stuff to start via slurm
eval "$(conda shell.bash hook)"
conda activate gene_phylogeny

# Main Program
file=$"WSC4_uniprot_e_1e-50_msa"

clustalo -i ../data/WSC4/${file}.afa \
	 --distmat-out ../data/WSC4/${file}_dist_matrix.csv \
         --percent-id \
	 -o /dev/null \
         --threads 36 \
         -v \
         --force \
	 --full \
         
