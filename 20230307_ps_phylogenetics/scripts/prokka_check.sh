#!/bin/bash
# Script to check that all files have a .gff file in them

genomes=$( find ../data/prokka -type d )

echo "$( echo "${genomes}" | wc -l) genomes found"

for genome in $genomes
do
	genome_name=${genome##*/}

	if [[ ! -f ${genome}/${genome_name}.gff ]]
	then
	
	echo "${genome_name} not translated"

	fi


done
