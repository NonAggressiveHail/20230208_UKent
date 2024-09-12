#!/bin/bash

#SBATCH --job-name=orthomcl_blast
#SBATCH --mail-type=END
#SBATCH --mail-user=jgh8@kent.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=14000
#SBATCH --error=./orthomcl_blast_logs/%a.out
#SBATCH --out=./orthomcl_blast_logs/%a.out
#SBATCH --time=48:00:00
#SBATCH --array=1-60%30

# conda needs special stuff to start via slurm
source /shared/home/jgh8/miniconda3/etc/profile.d/conda.sh

conda activate blast

# Only do setup on arrays 1 and 2
if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]
then

	
	echo 
	echo "$(date): Running setup on array $SLURM_ARRAY_TASK_ID"
	
	# make blast database if not already found
	if [[ ! -d "../data/orthomcl/goodProteins_db/" ]]
	then
	
		echo 
		echo "Running blastmakedb on $SLURM_ARRAY_TASK_ID"
		
		makeblastdb -in ../data/orthomcl/goodProteins.fasta \
					-out ../data/orthomcl/goodProteins_db/goodProteins \
					-dbtype prot
	else
		
		echo 
		echo "blastdb already found, skipping"
				
	fi

	# Make subfiles to be run in blast	
	echo 
	echo "Running seqkit split2 on array $SLURM_ARRAY_TASK_ID"
	
	seqkit split2 ../data/orthomcl/goodProteins.fasta \
				-p 60 \
				-f \
			   --out-dir ../data/orthomcl/split_sequences
			   
	echo "$(date): Done!"

else 

	echo 
	echo "$(date): Sleeping array $SLURM_ARRAY_TASK_ID whilst setup is done (5 minutes)"
	sleep 5m
	
fi

# Convert slurm array task id so it has leading zeros
file_number=$( printf '%03d' $SLURM_ARRAY_TASK_ID )

# Run blast
echo 
echo "$(date): Starting blastp on array $SLURM_ARRAY_TASK_ID"

blastp -query ../data/orthomcl/split_sequences/*${file_number}* \
	   -db ../data/orthomcl/goodProteins_db/goodProteins \
	   -out ../data/orthomcl/all_v_all_blast_${file_number}.tsv \
	   -outfmt 6 \
	   -num_threads ${SLURM_CPUS_PER_TASK}

echo	   
echo "$(date): Array $SLURM_ARRAY_TASK_ID complete, for orthomcl you will want to concatenate the blast outputs together"
echo  "then run orthomclBlastParser"
