#!bin/bash
#TODO add a project option for some namings and ANVIO
##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "Identifies the core genome in a selection of fasta files and"
	echo "produces a phylogenetic tree of it using anvi'o"
	echo 
	echo "options:"
	echo "  -h  display this help"
	echo "  -g  directory genomes are in"
	echo "  -p  pattern files have"
	echo "  -c  do genomes need decompressing?"
	echo "  -n  number of threads to use for parallel processing"
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
while getopts ":hg:p:cn:o:" option;
do
	case $option in
		h) # Display help
		  Help
		  exit;;
		g) # Directory genomes are in
		  gene_dir=$OPTARG;;
		p) # Pattern used to find files
	      file_pattern=$OPTARG;;
		c) # Do genomes need decompressing?
		  decompress=TRUE;;
		n) # Number of CPUs available
		  n_cpu=$OPTARG;;
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
echo "-c = ${decompress}"
echo "-n = ${n_cpu}"
echo "-o = ${out_dir}"
echo ""
echo "running anvio_phylogenomis on genomes in ${gene_dir},"
echo "using ${n_cpu} cpu's, and outputting results in ${out_dir}."
echo ""
if [ $decompress = TRUE ]
  then 
    echo "genomes will be decompressed before running"
fi

if [ $decompress = FALSE ]
  then 
    echo "genomes will not be decompressed before running"
fi

echo ""

# Make project name (probably want to redo my options for this...)
project=$(basename ${out_dir})
echo "project is ${project}"

# make output folders in data so I don't have to use raw_data
n_genomes=0
for file in $(find $gene_dir -type f -name $file_pattern)
do
   
  file_name=${file##*/}
  genome_name=${file_name%.fna*}
  
  mkdir -p ${out_dir}genomes/${genome_name}
 
  let n_genomes=n_genomes+1

done

echo "made $n_genomes folders"

## Reformat fasta files so anvio is happy with them
#if [ $decompress = TRUE ] 
#then
#  echo "Decompressing genomes"
#
#  find $gene_dir -type f -name $file_pattern |
#  parallel -j $n_cpu 'gzip -d {}' 
#
#  echo "Reformating files"
#  find $gene_dir -type f -name ${file_pattern%.gz} | 
#  eval "parallel -j $n_cpu 'anvi-script-reformat-fasta {} \
#                                                       -o ${out_dir}genomes/{/.}/{/.}_genomic.fna \
#                                                       --simplify-names' \
#													   --min-len 1000 \
#													   --seq-type NT"
#
#  echo "Recompressing genomes"
#
#  find $gene_dir -type f -name ${file_pattern%.gz} | 
#  parallel -j $n_cpu 'gzip {}'
#
#else
# 
#  echo "Reformating files"
#  find $gene_dir -type f -name $file_pattern | 
#  eval "parallel -j $n_cpu 'anvi-script-reformat-fasta {} \
#                                                       -o ${out_dir}genomes/{/.}/{/.}_genomic.fna \
#                                                       --simplify-names' \
#													   --min-len 1000 \
#													   --seq-type NT"
#fi
# 
#
## convert files into contigs.db files as this is the format anvio wants
#echo "Converting to contigs.db files"
#find ${out_dir}genomes/* -type d | 
#eval "parallel -j $n_cpu 'anvi-gen-contigs-database -f {}/{/}_genomic.fna \
#                                                    -o {}/{/}_genomic.db \
#                                                    -T 1 \
#											        -n ${project}'"
#
#echo "Done!"
#
## find hmms in contigs .db files
#echo "Finding hmms in contigs.db files"
#find ${out_dir}genomes/* -type d | 
#parallel -j $n_cpu 'anvi-run-hmms -c {}/{/}_genomic.db'
#echo "Done!"
#
## make file which says where the genomes are (anvi'o needs this later)
#echo "Making gene map file"
#echo -e "name\tcontigs_db_path" > ${out_dir}${project}_genome_storage.txt
#
#for dir in $(find ${out_dir}genomes -mindepth 1 -maxdepth 1 -type d)
#do
#  dirname=${dir##*/}
#  name=${dirname%%.*}
#  path=./genomes/${dirname}/${dirname}_genomic.db
#  line="${name}\\t${path}"
#
#  echo -e $line >> ${out_dir}${project}_genome_storage.txt
#
#done
#
#anvi-gen-genomes-storage -e ${out_dir}${project}_genome_storage.txt \
#			             -o ${out_dir}${project}-GENOMES.db
#						 
## above tested working, below not tested
## Run pangenome analysis 
anvi-pan-genome \
  -g ${out_dir}${project}-GENOMES.db \
  --output-dir ${out_dir} \
  --project-name ${project} \
  --num-threads ${n_cpu} \
  --minbit 0.5 \
  --mcl-inflation 10 \
  --min-occurrence 2


# Calculate ANI
#anvi-compute-genome-similarity -e ${out_dir}${project}_genome_storage.txt \
#                               -o ${out_dir}ANI \
#                               -p ${out_dir}${project}-PAN.db \
#     	                       -T ${n_cpu} \
#                               --program fastANI

# Get fasta sequencs for single copy core genes (SCGs) to later perform phylogenetics
#anvi-get-sequences-for-gene-clusters -p ${out_dir}${project}-PAN.db \
#                                     -g ${out_dir}${project}-GENOMES.db \
#                         		 --min-num-genomes-gene-cluster-occurs $(find ${out_dir}genomes -mindepth 1 -maxdepth 1 -type d | wc -l) \
#                                     --max-num-genes-from-each-genome 1 \
#                                     --concatenate-gene-clusters \
#                                     --output-file ${out_dir}${project}-SCGs.fna 

# Clean up fasta file with trimal
#trimal -in ${out_dir}${project}-SCGs.fna \
#       -out ${out_dir}${project}-SCGs_clean.fna \
#       -gt 0.5

# Generate phylogenetic tree with IQtree
#iqtree -s ${out_dir}${project}-SCGs_clean.fna \
#       -m WAG \
#       -bb 1000 \
#       -nt AUTO \
#       --threads-max ${n_cpu}

# generate layers order file, which organises the genomes in the pangeonome graph
#echo -e "item_name\tdata_type\tdata_value" > ${out_dir}${project}_phylogenomic_layer_order.txt

# add the newick tree as an order
#echo -e "SCGs_Bayesian_Tree\tnewick\t`cat ${out_dir}${project}-SCGs_clean.fna.contree`" \
#        >>  ${out_dir}${project}_phylogenomic_layer_order.txt
		
# import the layers order file
#anvi-import-misc-data -p ${out_dir}${project}-PAN.db \
#                      -t layer_orders ${out_dir}${project}_phylogenomic_layer_order.txt

echo ""
echo "Run Complete!"
 

