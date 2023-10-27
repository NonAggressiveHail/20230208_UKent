#!bin/bash

##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
 	echo "Identifies common hmms in genomes and "	
	echo "produces a phylogenetic tree of them using anvi'o"
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

if [ $decompress = TRUE ]
  then 
    echo "genomes will be decompressed before running"
fi

if [ $decompress = FALSE ]
  then 
    echo "genomes will not be decompressed before running"
fi

echo ""

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

# Reformat fasta files so anvio is happy with them
# This uses some reccomened anvio defaults
if [ $decompress = TRUE ] 
then
  echo "Decompressing genomes"

  find $gene_dir -type f -name $file_pattern |
  parallel -j $n_cpu 'gzip -d {}' 

  echo "Reformating files"
  find $gene_dir -type f -name ${file_pattern%.gz} | 
  eval "parallel -j $n_cpu 'anvi-script-reformat-fasta {} \
                                                       -o ${out_dir}genomes/{/.}/{/.}_genomic.fna \
                                                       --simplify-names' \
													   --min-len 1000 \
												   --seq-type NT"

  echo "Recompressing genomes"

  find $gene_dir -type f -name ${file_pattern%.gz} | 
  parallel -j $n_cpu 'gzip {}'

else
 
  echo "Reformating files"
  find $gene_dir -type f -name $file_pattern | 
  eval "parallel -j $n_cpu 'anvi-script-reformat-fasta {} \
                                                       -o ${out_dir}genomes/{/.}/{/.}_genomic.fna \
                                                       --simplify-names' \
													   --min-len 1000 \
												   --seq-type NT"
fi
 

# convert files into contigs.db files as this is the format anvio wants
echo "Converting to contigs.db files"
find ${out_dir}genomes/* -type d | 
parallel -j $n_cpu 'anvi-gen-contigs-database -f {}/{/}_genomic.fna \
                                              -o {}/{/}_genomic.db \
                                              -T 1'

echo "Done!"

# find hmms in contigs .db files
echo "Finding hmms in contigs.db files"
find ${out_dir}genomes/* -type d | 
parallel -j $n_cpu 'anvi-run-hmms -c {}/{/}_genomic.db'
echo "Done!"

# make file which says where the genomes are (anvi'o needs this later)
echo "Making gene map file"
echo -e "name\tcontigs_db_path" > ${out_dir}genome_storage.txt

for dir in $(find ${out_dir}genomes -mindepth 1 -maxdepth 1 -type d)
do
  dirname=${dir##*/}
  name=$(echo $dirname | cut -d _ -f1-2)
  path=./genomes/${dirname}/${dirname}_genomic.db
  line="${name}\\t${path}"

  echo -e $line >> ${out_dir}genome_storage.txt

done

# Get and concattenate the genes which are shared by all bacteria
echo "getting sequences for all ribosomal proteins in the Bacteria_71 gene set"
echo "this is what the anvio guys use on their website but may need optimising"
echo "for example, you may want to get SCGs if genomes are closely related"

anvi-get-sequences-for-hmm-hits --external-genomes ${out_dir}genome_storage.txt \
                                -o ${out_dir}concatenated_proteins.faa \
                                --hmm-source Bacteria_71 \
                                --gene-names ../raw_data/ribosomal_gene_names.txt \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate \

# Make phylogenetic tre from thesse hits
echo "Making phylogenetic tree..."

anvi-gen-phylogenomic-tree -f ${out_dir}concatenated_proteins.faa \
                           -o ${out_dir}phylogenomic_tree.tree

# Make an annotation file for all the genomes
echo "Making geneome tree annotations file"
echo -e "genome_id\tgenome" > ${out_dir}tree_annotations.txt

for dir in $(find ${out_dir}genomes -mindepth 1 -maxdepth 1 -type d)
do
 
  dirname=${dir##*/}
  name=$(echo $dirname | cut -d _ -f1-2)
  line="${name}\\t${name}"

  echo -e $line >> ${out_dir}tree_annotations.txt
done

echo ""
echo "Done!"
echo ""
echo "Tree can be viewed with anvi-interactive -p ${out_dir}phylogenomic-profile.db  --manual"
echo "                                         -t ${out_dir}phylogenomic_tree.tree"
echo "                                         -d ${out_dir}tree_annotations.txt"
echo "                                         --title mytitle"
echo "                                         --server-only"
echo "                                         --manual"

