#!/bin/bash
#TODO rewrite so that it just saves out all the seqkit grep commands, then I can run with slurm parallel
# ~ or just use GNU parallel for speedups on tesla? not sure 
#TODO add checks to skip already complete bits

##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "Given the blastp_annotations"
	echo "extracts loci ids for each siderophore protein per organism"
	echo "then extracts these loci with seqkit grep for use with orthomcl"
	echo 
	echo "options:"
	echo "  -j  blastp_annotation file" 
	echo "      produced by blastp_annotator.R"
	echo "  -p  directory with prokka .faa files in"
	echo "      this expects each .faa file to be in a folder with the organims name"
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
		j) # blast file
		  blast_file=$OPTARG;;
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
#   ~ extract loci sequences from faa file
# ~ done
#
# ~ compile together all the sequences from each organism, for each PVD number 

# Check all required args are present
if [ -z $blast_file ] || [ -z $prokka_dir ] || [ -z $out_dir ] || [ -z $threads ]
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

# append / to prokka_dir if not present
if [[ $prokka_dir != */ ]]
then
	prokka_dir="${prokka_dir}/"
fi

# Check args are what I expect
echo "blast_file: $blast_file"
echo "prokka_dir: $prokka_dir"
echo "out_dir: $out_dir"
echo "threads: $threads"
echo



# Run rscript to split out all my loci into files
echo 
echo "Running loci_collector.R"
echo "Command: Rscript loci_collector.R ${blast_file} ${out_dir} ${threads}"

#Rscript loci_collector.R ${blast_file} ${out_dir} ${threads}

## Run seqkit
echo
echo "Running seqkit"
echo


n_organisms=$(ls ${out_dir} | wc -l)
counter=0

# Use seqkit to pull out the loci sequences
for dir in ${out_dir}*
do
	#echo "dir is $dir"
	organism=${dir##*/}
	#echo "organism is $organism"
	
	for sid_loci_file in ${out_dir}${organism}/*_loci_ids.txt
	do
		# TODO change this so I dont have "_loci_ids_sequences" in the sid_id	
		#echo  "sid_loci_file is $sid_loci_file"
		sid_id=$(basename $sid_loci_file .txt)
		#echo "sid_id is $sid_id"
		
		#echo "Seqkit command:"
		#echo 
		#echo "seqkit grep \ "
		#echo "-f $sid_loci_file \ "
		#echo "-o ${out_dir}${organism}/${sid_id}_sequences.faa \ "
		#echo "-j $threads \ "
		#echo "--quiet \ "
		#echo "${prokka_dir}${organism}/${organism}.faa"		
		
		seqkit grep \
			-f $sid_loci_file \
			-o ${out_dir}${organism}/${sid_id}_sequences.faa \
			-j $threads \
			--quiet \
			${prokka_dir}${organism}/${organism}.faa
			
		
	done
	
	# Progress counter
	let counter=counter+1
	per_comp=$(echo "scale=2 ; $counter * 100 / $n_organisms" | bc)
	echo "${per_comp}% complete"
	
done


# Gather all the seqeuenees together in the same order for each
echo 
echo "Gathering all sequences together"
echo

all_sid_ids=$( tail -n +2 ${blast_file} |
			   xsv select 2 |
			   sort |
			   uniq )

counter=0
for dir in ${out_dir}*
do

	organism=${dir##*/}
	echo "Working on ${organism}"

	# Make output file
	echo "Make output file"
	echo "> ${organism}_siderophore_genes" > ${dir}/concatenated_sid_genes.faa

	echo "copying sid genes"
	echo

	for sid_id in $all_sid_ids
	do
		#echo "grep command is grep -v '>' ${dir}/${sid_id}_loci_ids_sequences.faa >> ${dir}/concatenated_sid_genes.faa"

		# this should append but also remove trailing newlines and whitespace
		echo -n `grep -v ">" ${dir}/${sid_id}_loci_ids_sequences.faa` | tr -d "[:space:]" >> ${dir}/concatenated_sid_genes.faa

	done

done

echo "Script complete"

