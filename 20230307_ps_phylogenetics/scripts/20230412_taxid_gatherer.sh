#!/bin/bash
# TODO change to using epost, it was erroring before as I had "OX" at the end of the input file for some reason
#######################################################
# help
#######################################################

Help()
{
  # Display Help
  echo ""
  echo "This takes a list of NCBI taxonomy UID numbers (TXID), and gathers"
  echo "the full taxonomy for each"
  echo
  echo "Syntax: 20230303_protein_uid_to_taxonomy.sh [-h|-i]"
  echo "options:"
  echo "  -h  Print this Help."
  echo "  -i  Input NCBI UID list"
  echo "  -o  Output file location"
  echo
  echo "Author: Jake Hudson"
  echo
}

#####################################################
# main program 
#####################################################
#####################################################

# process input optios
# Get the options
while getopts ":hi:o:" option;
do
  case $option in
    h) # display Help
      Help
      exit;;
    i) # Input 
      in=$OPTARG;; 
    o) # Output
      out=$OPTARG;;
   \?) #Invalid options
      echo "Error: Invalid option"
      exit;;
  esac
done

###################################################
# Main Program
###################################################
echo "taxid,name,kingdom,phylum,class,order,family,genus" > ${out}

for id in $(cat $in)
  do
    echo "Searching for ${id}"

    efetch -db taxonomy -id ${id} -format xml |
      xtract -pattern Taxon -tab "," -first TaxId ScientificName \
      -group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -GNUS "(-)" \
      -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName \
      -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName \
      -block "*/Taxon" -match "Rank:class" -CLSS ScientificName \
      -block "*/Taxon" -match "Rank:order" -ORDR ScientificName \
      -block "*/Taxon" -match "Rank:family" -FMLY ScientificName \
      -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName \
      -group Taxon -tab "," -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&GNUS" >> ${out}
       
       sleep 0.3s
done


