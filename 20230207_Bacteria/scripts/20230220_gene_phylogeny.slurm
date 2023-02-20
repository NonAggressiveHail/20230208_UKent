#!/bin/bash

#SBATCH --job-name=WSC4_phylogeny
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --partition=icarus-biol-cpu-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36

20230209_gene_phylogeny.sh -f ../data/WSC4/
