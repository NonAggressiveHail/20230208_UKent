#!/bin/bash

#######################################################
# help
#######################################################

Help()
{
  # Display Help
  echo "Add description of the script functions here."
  echo
  echo "Syntax: scriptTemplate [-g|h|v|V]"
  echo "options:"
  echo "g     Print the GPL license notification."
  echo "h     Print this Help."
  echo "v     Verbose mode."
  echo "V     Print software version and exit."
  echo
}

#####################################################
# main program 
#####################################################
#####################################################

# process input optios
# Get the options
while getopts ":hf:" option;
do
  case $option in
    h) # display Help
      Help
      exit;;
    f) # Enter a file to run on
       file=$OPTARG;; 
   \?) #Invalid options
      echo "Error: Invalid option"
      exit;;
  esac
done

###################################################
# Main Program
###################################################
#echo genome
echo "I am running on ${file}"

#get accession number
tmp=${file%_genomic.*}
accession=${tmp##*/}
echo "accession is ${accession}"

#get directory to save outputs too
dir=${file%/*}

#run prokka
prokka --outdir ${dir}prokka --prefix ${accession} ${file} 






