#!/bin/bash
# Check how many characters are in a filename
files=$( find ../raw_data/genomes -type f -name *.fna )

for file in $files
do

	filename=${file##*/}
	file_title=${filename%.fna}

	n_characters=$(echo ${file_title} | wc -c)
	
	if [[ $n_characters -gt 30 ]]
	then
	
	echo "${file_title} has ${n_characters} characters"
	
	
	fi

done
