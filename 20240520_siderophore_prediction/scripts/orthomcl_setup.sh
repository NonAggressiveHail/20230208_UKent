#!/bin/bash

##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "Given a folder with prokka results in, this script will prepare the fasta files for use with orthomcl"
	echo 
	echo "options:"
	echo "  -p  directory with prokka .faa files in"
	echo "      this expects each .faa file to be in a folder with the organims name"
	echo "  -o  output directory"
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
while getopts ":p:o:h" option;
do
	case $option in
		p) # prokka directory
		  prokka_dir=$OPTARG;;
		o) # output directory
		  out_dir=$OPTARG;;
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
# Check inputs
echo "Checking inputs"
echo "prokka_dir: ${prokka_dir}"
echo "out_dir: ${out_dir}"

if [ -z ${prokka_dir} ] || [ -z ${out_dir} ]
then
    echo "Error: Missing required inputs"
    Help
    exit
fi

# Remove trailing slashes
prokka_dir=$(echo ${prokka_dir} | sed 's:/*$::')
out_dir=$(echo ${out_dir} | sed 's:/*$::')

# Set so that script exits on errir
set -e

# Setup for Orthomcl
## Ensure fasta files are compliant
### Genate file with a 3 digit code for each bacteria

# Make output location
mkdir -p ${out_dir}/compliantFasta/

## Make a temporary dir to store the compliant fasta files
mkdir -p ./temp/
cd ./temp/

echo -e "name,code" > ../${out_dir}/taxon_codes.csv
organism_number=1

# Run orthomclAdjustFasta
n_organisms=$(ls ../${prokka_dir} | wc -l)
for dir in ../${prokka_dir}/*
do 

	# Make taxon code/organism translation table
	organism=${dir##*/}
	taxon_code=$(printf "%03d" ${organism_number})
	

	echo -e "${organism},${taxon_code}" >> ../${out_dir}/taxon_codes.csv
	
	# Run orthomclAdjustFasta
	orthomclAdjustFasta ${taxon_code} ${dir}/proteins.faa 1 
		
	per_comp=$(echo "scale=2 ; $organism_number * 100 / $n_organisms" | bc)
	echo "${per_comp}% complete"
	
	let organism_number=organism_number+1
	
done

# Move all compliant fastas to outdir, as we can't control the output of orthomclAdjsutFAsta
mv *.fasta ../${out_dir}/compliantFasta/
cd ../
rmdir ./temp/

# Run orthomcl filterfasta
cd ${out_dir}
orthomclFilterFasta ./compliantFasta 10 20 
cd -

echo 

