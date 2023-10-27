#!/bin/bash

#SBATCH --job-name=p_aeruginosa_phylogeny
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=biosoc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=14000
#SBATCH --error=p_aeruginosa_phylogeny.err
#SBATCH --out=p_aeruginosa_phylogeny.out
#SBATCH --time=24:00:00

# conda needs special stuff to start via slurm
eval "$(/shared/home/jgh8/miniconda3/bin/conda shell.bash hook)"

conda activate anvio-7.1

# Main Program
bash 20230307_anvio_phylogenomics.sh -g ../raw_data/p_aeruginosa_genomes/ \
                                     -p "*.fna.gz" \
                                     -n $SLURM_CPUS_PER_TASK \
                                     -o ../data/p_aeruginosa_phylogeny/ \
                                     -c 
