#!/bin/bash
#TODO rewrite so that it just saves out all the seqkit grep commands, then I can run with slurm parallel
##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "Given the compiled_jackhmmer_data_corrected_names.csv"
	echo "extracts loci ids for each siderophore protein per organism"
	echo "then extracts these loci with seqkit grep for use with orthomcl"
	echo 
	echo "options:"
	echo "  -j  jackhmmer file" 
	echo "      produced by jackhmmer_compiler.R"
	echo "  -p  directory with prokka .ffn files in"
	echo "      this expects each .ffn file to be in a folder with the organims name"
	echo "  -o  output directory"
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
while getopts ":j:p:o:t:h" option;
do
	case $option in
		j) # jackhmmer file
		  jackhmmer_file=$OPTARG;;
		p) # prokka directory
		  prokka_dir=$OPTARG;;
		o) # output directory
		  out_dir=$OPTARG;;
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
# ~ filter compiled_jackhmmer_data_corrected_names.csv to only have PVD enzymes
#
# ~ for each target organism in compiled_jackhmmer_data_corrected_names.csv
#   ~ get loci numbers for each PVD enzyme
#   ~ extract loci sequences from fna file
# ~ done
#
# ~ compile together all the sequences from each organism, for each PVD number 

# Run rscript to split out all my loci into files
echo 
echo "Running loci_collector.R"
echo 

Rscript loci_collector.R ${jackhmmer_file} ${out_dir}/fasta_files/ ${threads}

# Run seqkit
echo
echo "Running seqkit"
echo


n_organisms=$(ls ${out_dir}/fasta_files/ | wc -l)
counter=0
for dir in ${out_dir}/fasta_files/*
do
	#echo "dir is $dir"
	organism=${dir##*/}
	#echo "organism is $organism"
	
	for sid_loci_file in ${out_dir}/fasta_files/${organism}/*_loci_ids.txt
	do
		#echo "sid_loci_file is $sid_loci_file"
		sid_id=$(basename $sid_loci_file .txt)
		#echo "sid_id is $sid_id"
		
		#echo "Seqkit command:"
		#echo 
		#echo "seqkit grep \ "
		#echo "-f $sid_loci_file \ "
		#echo "-o ${out_dir}/fasta_files/${organism}/${sid_id}_sequences.fna "
		#echo "   ${prokka_dir}/${organism}/${organism}.ffn"		
		
		seqkit grep \
			-f $sid_loci_file \
			-o ${out_dir}/fasta_files/${organism}/${sid_id}_sequences.fna \
			-j $threads \
			--quiet \
			${prokka_dir}/${organism}/${organism}.ffn
			
		
	done
	
	# Progress counter
	let counter=counter+1
	per_comp=$(echo "scale=2 ; $counter * 100 / $n_organisms" | bc)
	echo "${per_comp}% complete"

done

echo "Script complete"
