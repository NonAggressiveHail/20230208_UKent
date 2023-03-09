#!bin/bash

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
	echo "  -c  do genomes need decompressing?"
	echo "  -n  number of threads to use for parallel processing"
	echo "  -o directory to output results in"
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
while getopts ":hg:c:n:o:" option;
do
	case $option in
		h) # Display help
		  Help
		  exit;;
		g) # Directory genomes are in
		  gene_dir=$OPTARG;;
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
# make output folders in data so I don't have to use raw_data
for dir in $(find $gene_dir -mindepth 1 -maxdepth 1 -type d)
do

  dirname=${dir##*/}
  mkdir -p ${out_dir}genomes/${dirname}

done


# decompress files if needed
if [ decompress = TRUE ] 
  then
    
    echo "Decompressing"

    find $gene_dir -mindepth 1 -maxdepth 1 -type d |
      parallel -j $n_cpu 'gzip -d {}/{/}_genomic.fna.gz'

    echo "Done!"

fi


# Reformat fasta files so anvio is happy with them
echo "Reformatting files"
find $gene_dir -mindepth 1 -maxdepth 1 -type d | 
 eval "parallel -j $n_cpu 'anvi-script-reformat-fasta {}/{/}_genomic.fna \
                                                 -o ${out_dir}genomes/{/}/{/}_genomic.fna \
                                                 --simplify-names'"

# convert files into contigs.db files as this is the format anvio wants
echo "Converting to contigs.db files"
find $gene_dir -mindepth 1 -maxdepth 1 -type d | 
 eval "parallel -j $n_cpu 'anvi-gen-contigs-database -f ${out_dir}genomes/{/}/{/}_genomic.fna \
                                                     -o ${out_dir}genomes/{/}/{/}_genomic.db \
                                                     -T 1'"

echo "Done!"

# find hmms in contigs .db files
# I'm not sure why we do this but its in the tutorial
echo "Finding hmms in contigs.db files"
find $gene_dir -mindepth 1 -maxdepth 1 -type d | 
 eval "parallel -j $n_cpu 'anvi-run-hmms -c ${out_dir}genomes/{/}/{/}_genomic.db'"
echo "Done!"

# make file which says where the genomes are (anvi'o needs this later)
echo "Making gene map file"
echo -e "name\tcontigs_db_path" > ${out_dir}genome_storage.txt

for dir in $(find ${out_dir}genomes -mindepth 1 -maxdepth 1 -type d)
do

  dirname=${dir##*/}
  name=${dirname%%.*}
  path=./genomes/${dirname}/${dirname}_genomic.db
  line="${name}\\t${path}"

  echo -e $line >> ${out_dir}genome_storage.txt

done

# Get and concattenate the genes which are shared by all bacteria
echo "Getting sequences for all shared hmms, you may want to change this depending on you question"

anvi-get-sequences-for-hmm-hits --external-genomes ${out_dir}genome_storage.txt \
                                -o ${out_dir}concatenated_proteins.faa \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate \

# Make phylogenetic tre from thesse hits
echo "Making phylogenetic tree..."

anvi-gen-phylogenomic-tree -f ${out_dir}concatenated_proteins.faa \
                           -o ${out_dir}phylogenomic_tree.tree

echo "Done!"

 

