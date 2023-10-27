#!bin/bash
##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "given a folder, moves each file in it" 
	echo "into its own subdirectory" 
	echo "options:"
	echo "  -i input files"
	echo "  -t threads to use"
	echo "  -h display this help"
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
while getopts ":hi:t:" option;
do
	case $option in
		h) # Display help
			Help
			exit;;
		i) # inputs
		  dir=$OPTARG;;
		t) # threads
		  threads=$OPTARG;;		  
		\?) #Invalid option
			echo ""
			echo "Error: Invalid option"
			Help
			exit;;
		esac
done

###############################
# Main Program
###############################
# Debug stuff
echo ""
echo "-i ${dir}"
echo "-t ${threads}"
echo ""

# Count files I have and confirm if we want to proceede 
count=$(find $dir -maxdepth 1 -type f | wc -l)

while true
	do
		read -p "${count} files found, continue? y/n: " yn
		case $yn in 
			[Yy]* ) # Proceede with download
			          echo "Making folders now"
			          break;;
			[Nn]* ) # Abort Download
			          echo "Aborting"
			          exit;;
			* ) #Invalid options
			    echo "Please answer y or n"
		esac
done

# make folders and move folders in parallel
find $dir -maxdepth 1 -type f | 
  parallel -j $threads --plus 'mkdir {%_genomic.fna.gz} 
                               mv {} {%_genomic.fna.gz}/{/}'

















