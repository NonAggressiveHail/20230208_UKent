#!/bin/bash

#SBATCH --job-name=WSC4_phylogeny
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --error=WSC4_phylogeny.err
#SBATCH --out=WSC4_phylogeny.out
#SBATCH --time=03:00:00

#conda needs special stuff to start via slurm
eval "$(conda shell.bash hook)"
conda activate gene_phylogeny

bash 20230209_gene_phylogeny.sh -f ../raw_data/proteins/WSC4.fna -b ../raw_data/uniprot/uniprot_all.fasta -o ../data/WSC4/WSC4_uniprot -t -c 36 
