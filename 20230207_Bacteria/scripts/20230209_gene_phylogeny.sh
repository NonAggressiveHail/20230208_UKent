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
	echo 
}

#############################################
# Process input options
#############################################

# Get options
while getopts ":hf:" option;
do
	case $option in
		h) # display Help
			Help
			exit;;
		f) # file to be run on
			file=$OPTARG;;
		\?) # Invalid option
			echo "Error: Invalid option"
			Help
			exit;;
		esac
done

##############################################################
# Main Program
#############################################################



