#!/bin/bash


##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "Given fasta files pulls out sequences containing a specific string"
	echo "and rotates them to have a common start sequence"
	echo 
	echo "options:"
	echo "  -f  directory with fna files" 
	echo "  -o  output directory"
    echo "  -n  nucleotide sequence to start from"
    echo "  -s  pattern in headers"
	echo "  -t  number of threads"
	echo "  -h  display this help"
	echo
	echo "Author: Jake Hudson"
	echo
}

#################################
# Setup
#################################


#################################
# Process input options
#################################

# Get Options
while getopts ":f:o:n:s:t:h" option;
do
	case $option in
		f) # file containing genomes
		  genomes_dir=$OPTARG;;
		o) # output directory
		  out_dir=$OPTARG;;
        n) # nucleotide sequence to start from
          nucleotide=$OPTARG;;
		s) # search string
		  search_string=$OPTARG;;
		t) # number of threads
		  threads=$OPTARG;;
		h) # Help
		  Help
		  exit;;
		\?) #Invalid option
		  echo "Error: Invalid option"
		  Help
		  exit;;
		esac
done

###############################
# Main Program
###############################
# steps
# ~ pull out seqnencs with headers matching string into one file
# ~ rotate sequences to have a common start sequence
# ~ save into seperate files per sequence with other sequences from the original file

# Check all required args are present
if [ -z $genomes_dir ] || [ -z $out_dir ] || [ -z $search_string ] || [ -z $threads ] || [ -z $nucleotide ]
then
    echo "Error: Missing required arguments"
    Help
    exit
fi

# append / to out_dir if not present
if [[ $out_dir != */ ]]
then
	out_dir="${out_dir}/"
fi

# append / to genomes_dir if not present
if [[ $genomes_dir != */ ]]
then
	genomes_dir="${genomes_dir}/"
fi

# Echo arguments
echo "Genomes directory: $genomes_dir"
echo "Output directory:  $out_dir"
echo "Search string:     $search_string"
echo "Threads:           $threads"
echo "Nucleotide:        $nucleotide"

# Make output directory
mkdir -p $out_dir

# Make table of filenames with acccession numbers

#echo "filename,accession" > $out_dir"filenames_accessions.txt"
#for( file in $genomes_dir*.fna )
#do
#    filename=$(basename $file)
#    accession=$(grep -oP '^>[^|]+' $file) # this may match multiple lines, try and then check
#
#    echo "$filename,$accession" >> $out_dir"filenames_accessions.txt"
#done

# Combine sequences together
echo $search_string > ids.txt
seqkit grep ${genomes_dir}*.fna -n -f ids.txt -r -o ${out_dir}unaligned_sequences.fna
rm ids.txt

# Rotate sequences together
rotate  -s $nucleotide -o ${out_dir}rotated_sequences.fna ${out_dir}unaligned_sequences.fna

# recombine rotated sequences with original sequences
# pull out lines which dont match chromomosome one at at time?