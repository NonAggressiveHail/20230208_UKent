#!bin/bash
##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "Produces summary statistics of genomes using quast"
	echo 
	echo "options:"
	echo "  -h  display this help"
	echo "  -g  directory genomes are in"
	echo "  -p  pattern files have"
	echo "  -n  number of threads for each job"
	echo "  -j  number of parallel jobs"
	echo "  -o  directory to output results in"
	echo
	echo "Author: Jake Hudson"
	echo
}

#################################
# Setup
#################################
n_cpu=1
decompress=FALSE

#################################
# Process input options
#################################

# Get Options
while getopts ":hg:p:n:j:o:" option;
do
	case $option in
		h) # Display help
		  Help
		  exit;;
		g) # Directory genomes are in
		  gene_dir=$OPTARG;;
		p) # Pattern used to find files
	          file_pattern=$OPTARG;;
		n) # Number of CPUs/job
		  n_cpu=$OPTARG;;
		j) # Number of jobs
		  n_job=$OPTARG;;
		o) # output location
		  out_dir=$OPTARG;;
		\?) #Invalid option
			echo "Error: Invalid option"
			Help
			exit;;
		esac
done

###############################
# Main Program
###############################
# Some debugging info
echo ""
echo "-g = ${gene_dir}"
echo "-p = ${file_pattern}"
echo "-n = ${n_cpu}"
echo "-o = ${out_dir}"
echo ""
echo "running anvio_phylogenomis on genomes in ${gene_dir},"
echo "using ${n_cpu} cpu's, and outputting results in ${out_dir}."
echo ""

# find files
files=$(find $gene_dir -type f -name $file_pattern)

echo "Running on:"
echo "${files}"
echo ""

# make output folders in data so I don't have to use raw_data
#n_genomes=0
#for file in files
#do
#   
#  file_name=${file##*/}
#  genome_name=${file_name%.fna*}
#  
#  mkdir -p ${out_dir}genomes/${genome_name}
# 
#  let n_genomes=n_genomes+1
#
#done
#
#echo "made $n_genomes folders"

# Make output folder 
mkdir -p ${out_dir}QUAST

# Run quast
#quast.py \
#	-o ${out_dir}QUAST \
#	--threads $(( $n_cpu * $n_job )) \
#	--silent \
#	${files}

# Run prokka
# center and locus tag options are a bugfix
mkdir -p ${out_dir}prokka

#eval "parallel -j ${n_job} 'prokka --outdir ${out_dir}prokka/{/.} \
#	  		           --prefix {/.} \
#				   --genus Pseudomonas \
#	                           --species aeruginosa \
#         	                   --cpus ${n_cpu} \
#				   --compliant \
#			           --force \
#			           --centre C \
#				   --locustag L \
#	          	           {} '" ::: ${files}

# Run phigaro
## NB this uses an older version of pythin so needs its own conda
mkdir -p ${out_dir}phigaro
mkdir -p ${out_dir}temp_files

## First we need to remove contigs <= 20000 for phigaro to work
#eval "parallel -j ${n_job} 'reformat.sh in={} out=${out_dir}temp_files/{/.}_trim.fna -minlength=20000'" ::: ${files}

#eval "parallel -j ${n_job} 'phigaro -f {} \
#				    -e tsv \
#                                   -o ${out_dir}phigaro/{/.}_phages \
#				    -t ${n_cpu} \
#                                 --not-open '" ::: ${out_dir}temp_files/*trim.fna

#rm ${out_dir}temp_files/*_trim.fna

# Run platon
mkdir -p ${out_dir}platon

#eval "parallel -j ${n_job} 'platon --db ../programs/platon_db/ \
#				   --output ${out_dir}platon/{/.} \
#				   --threads ${n_cpu} \
#				   {} '" ::: ${files}  

# Run isescan
mkdir -p ${out_dir}isescan 

eval "parallel -j ${n_job} 'isescan.py --seqfile {} \
				       --output ${out_dir}isescan/{/.} \
				       --nthread ${n_cpu}'" ::: ${files}
























