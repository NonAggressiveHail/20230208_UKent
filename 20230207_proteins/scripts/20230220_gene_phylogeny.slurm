#!/bin/bash

#SBATCH --job-name=WSC4_phylogeny_e_1e-50
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=14G
#SBATCH --error=WSC4_phylogeny_e_1e-50.err
#SBATCH --out=WSC4_phylogeny_e_1e-50.out
#SBATCH --time=1:00:00

# conda needs special stuff to start via slurm
eval "$(conda shell.bash hook)"
conda activate gene_phylogeny

# Main Program
bash 20230209_gene_phylogeny.sh -f ../raw_data/proteins/WSC4.fna -b ../raw_data/uniprot/uniprot_all.fasta -o ../data/WSC4/WSC4_uniprot_e_1e-50 -t -c 36 -e 1e-50
