#!/bin/bash

#TODO change to using FastTree MP for multithreading
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
	echo "  -b    multifasta to search -f against"
	echo "  -o    outputs filename"
	echo "  -t    should a phylogenetic tree be built?"
	echo "  -c    number of cpu's for jackhmmer"
	echo "  -e    e-value for jackhmmer cutoff"
	echo 
}

#############################################
# Process input options
#############################################
# set variables
# fasttree off by default
tree=FALSE
# default jackhmmer e-value
e_value=10

# Get options
while getopts ":hf:b:o:tn:c:e:" option;
do
	case $option in
		h) # display Help
			Help
			exit;;
		f) # file to be run on
			file=$OPTARG;;
		b) #multifasta to search in
		  db=$OPTARG;;
		o) #output file location
		  out=$OPTARG;;
		t) # should fasttree be made from msa?
		   tree=TRUE;;
		c) # number of cpus to use for jackhmmer
		   cpu=$OPTARG;;
		e) # e-value for jackhmmer
		   e_value=$OPTARG;;
		\?) # Invalid option
			echo "Error: Invalid option"
			Help
			exit;;
		esac
done

##############################################################
# Main Program
#############################################################

#run jackhmmer
echo "running jackhmmer"
jackhmmer -o ${out}.txt -A ${out}_msa.sto -E $e_value --cpu ${cpu} ${file} ${db}

#run fasttree
if [ "$tree" = TRUE ] 
	then

	echo "converting ${out}_msa.sto to ${out}_msa.afa"
  	esl-reformat --informat stockholm -o ${out}_msa.afa afa ${out}_msa.sto

	echo "building tree from jackhmmer msa"
	#fasttree ${out}_msa.afa > ${out}.tree
	FastTreeMP -fastest ${out}_msa.afa > ${out}.tree

fi

# produce annotation files for iTOL
echo "Making iTOL annotations"
Rscript 20230216_iTOL_annotator.R ${out}_msa.afa ${out}

echo "Script complete!"





