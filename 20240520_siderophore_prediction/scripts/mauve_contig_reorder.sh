#!bin/bash

genomes=$( find ../raw_data/IATS_sequences/assemblies/ -type f -name *.fna )

for file in $genomes
do
	# Make output directories
	organism=$( basename $file .fna )
	output=../data/oriented_genomes/$organism 
	
	if [ -f $output/${organism}_reoriented.fasta ]
	then
	
		echo "${organism}_reoriented.fasta found" 
		
	else 
	
		echo "${organism}_reoriented.fasta not found, running Mauve" 
		
		rm -rf ../data/oriented_genomes/$organism
		mkdir -p ../data/oriented_genomes/$organism
	
		# Run Mauve contig aligner
		java -Xmx500m \
			-cp ../programs/mauve_snapshot_2015-02-13/Mauve.jar \
			org.gel.mauve.contigs.ContigOrderer \
			-output $output \
			-ref ../data/oriented_genomes/Pa_PAO1_107/Pa_PAO1_107_reoriented.fasta \
			-draft $file
		
		# Move final alignment to one file 
		newest=$(ls $output -1 -t | head -n 1)
		echo "move command:"
		echo "cp $output/$newest/${organism}.fna.fas $output/${organism}_reoriented.fasta"
		cp $output/$newest/${organism}.fna.fas $output/${organism}_reoriented.fasta
		
	fi 
	
	
	
done

echo "Complete!"



