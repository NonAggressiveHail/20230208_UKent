#!/bin/bash

###############################################
# Help
###############################################

Help()
{
	# Display Help
	echo "Given a gene, this script will search for related genes with
jackhmmmer, then produce a phylogentic tree of related genes"
	echo 
	echo "options:"
	echo "  -f    file to run on"
	echo "  -h    display this help"
	echo "  -b    directory to build database from"
	echo 
}

#############################################
# Process input options
#############################################

# Get options
while getopts ":hf:b:" option;
do
	case $option in
		h) # display Help
			Help
			exit;;
		f) # file to be run on
			file=$OPTARG;;
		b) #directory to build hmm database from
		  dir=$OPTARG;;
		\?) # Invalid option
			echo "Error: Invalid option"
			Help
			exit;;
		esac
done

##############################################################
# Main Program
#############################################################

#find cds files to unzip
echo "finding and unzipping files in ${dir}"
find $dir -type f -name *_cds_*.fna.gz -exec gzip -d -f -k {} +

#concatenate them all together so jackhmmer can search something
echo "finding file in ${dir} to build database from"
find $dir -type f -name *_cds_*.fna -exec cat {} + > ${dir}/database.fna

#run jackhmmer
echo "running jackhmmer"
jackhmmer -o ${dir}/jackhmmer.out -A ${dir}/jackhmmer_msa.fna $file ${dir}/database.fna






