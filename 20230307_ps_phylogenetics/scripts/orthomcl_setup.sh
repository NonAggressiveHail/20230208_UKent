#!/bin/bash
#TODO rewrite so that it just saves out all the seqkit grep commands, then I can run with slurm parallel
# ~ or just use GNU parallel for speedups on tesla? not sure 
#TODO add checks to skip already complete bits
#TODO check to see if jackhmmer can just output matches to file...
#TODO split into "Orthomcl prep" and "orthomcl run" so I can do an array submission

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
#   ~ extract loci sequences from faa file
# ~ done
#
# ~ compile together all the sequences from each organism, for each PVD number 

echo

# Run rscript to split out all my loci into files
#echo 
#echo "Running loci_collector.R"
#echo 
#
#Rscript loci_collector.R ${jackhmmer_file} ${out_dir}/fasta_files/ ${threads}
#
## Run seqkit
#echo
#echo "Running seqkit"
#echo
#
#
#n_organisms=$(ls ${out_dir}/fasta_files/ | wc -l)
#counter=0
#for dir in ${out_dir}/fasta_files/*
#do
#	#echo "dir is $dir"
#	organism=${dir##*/}
#	#echo "organism is $organism"
#	
#	for sid_loci_file in ${out_dir}/fasta_files/${organism}/*_loci_ids.txt
#	do
#		# TODO change this so I dont have "_loci_ids_sequences" in the sid_id	
#		#echo  "sid_loci_file is $sid_loci_file"
#		sid_id=$(basename $sid_loci_file .txt)
#		#echo "sid_id is $sid_id"
#		
#		#echo "Seqkit command:"
#		#echo 
#		#echo "seqkit grep \ "
#		#echo "-f $sid_loci_file \ "
#		#echo "-o ${out_dir}/fasta_files/${organism}/${sid_id}_sequences.faa "
#		#echo "   ${prokka_dir}/${organism}/${organism}.faa"		
#		
#		seqkit grep \
#			-f $sid_loci_file \
#			-o ${out_dir}/fasta_files/${organism}/${sid_id}_sequences.faa \
#			-j $threads \
#			--quiet \
#			${prokka_dir}/${organism}/${organism}.faa
#			
#		
#	done
#	
#	# Progress counter
#	let counter=counter+1
#	per_comp=$(echo "scale=2 ; $counter * 100 / $n_organisms" | bc)
#	echo "${per_comp}% complete"
#
#done

# Gather all the seqeuenees together
## Get Id of each siderophore protein
#all_sid_ids=$( xsv search -i -s 3 -n "Siderophore" $jackhmmer_file |
#               xsv select 6 |
#			   sort |
#			   uniq )
#
### Gather sequences into one file
#mkdir -p ${out_dir}/collated_sequences/
#
#for sid_id in $all_sid_ids
#do
#	
#	# delete file first to make sure we dont duplicate sequences
#	rm -f ${out_dir}/collated_sequences/${sid_id}_matches.faa
#	
#	find ${out_dir}/fasta_files/ -type f -name ${sid_id}_loci_ids_sequences.faa -exec cat {} >> ${out_dir}/collated_sequences/${sid_id}_matches.faa \;
#	
#	exit	
#	
#done
#
#echo "Script complete"

# Run Orthomcl
## Ensure fasta files are compliant
### Genate file with a 3 digit code for each bacteria

# Make output location
#mkdir -p ${out_dir}/compliantFasta/
#
## Make a temporary dir to store the compliant fasta files
#mkdir -p ./temp/
#cd ./temp/
#
#echo -e "name,code" > ../${out_dir}/taxon_codes.csv
#organism_number=1
#
## Run orthomclAdjustFasta
#n_organisms=$(ls ../${prokka_dir} | wc -l)
#for dir in ../${prokka_dir}/*
#do 
#
#	# Make taxon code/organism translation table
#	organism=${dir##*/}
#	taxon_code=$(printf "%03d" ${organism_number})
#	
#	echo -e "${organism},${taxon_code}" >> ../${outdir}/taxon_codes.csv
#	
#	# Run orthomclAdjustFasta
#	orthomclAdjustFasta ${taxon_code} ${dir}/${organism}.faa 1 
#		
#	per_comp=$(echo "scale=2 ; $organism_number * 100 / $n_organisms" | bc)
#	echo "${per_comp}% complete"
#	
#	let organism_number=organism_number+1
#		
#done
## Move all compliant fastas to outdir, as we can't control the output of orthomclAdjsutFAsta
#mv *.fasta ../${out_dir}/compliantFasta/
#cd ../
#rmdir ./temp/

## Run orthomcl filterfasta
cd ${out_dir}
orthomclFilterFasta ./compliantFasta 10 20 
cd -



