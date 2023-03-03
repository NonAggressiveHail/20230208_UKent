#!/bin/bash

#SBATCH --job-name=WSC4_phylogeny_e_0.01
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --error=WSC4_phylogeny_e_0.01.err
#SBATCH --out=WSC4_phylogeny_e_0.01.out
#SBATCH --time=24:00:00

# conda needs special stuff to start via slurm
eval "$(conda shell.bash hook)"
conda activate gene_phylogeny

# Main Program
bash 20230209_gene_phylogeny.sh -f ../raw_data/proteins/WSC4.fna -b ../raw_data/uniprot/uniprot_all.fasta -o ../data/WSC4/WSC4_uniprot_e_0.01 -t -c 36 -e 0.01
