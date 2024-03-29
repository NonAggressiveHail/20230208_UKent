#!bin/bash
# TODO change so that assembly report is also downloaded
# TODO make it so that each genome is in its own folder, so asseblb reorts can go there too
##################################
# Help()
##################################

Help()
{
	# Display Help
	echo 
	echo "Searches NCBI database for sequences and downloads them"
	echo 
	echo "options:"
	echo "  -s  search term." 
	echo "      you probably want foo [ORGN] AND latest [filter]"
	echo "      help here (https://linsalrob.github.io/ComputationalGenomicsManual/Databases/NCBI_Edirect.html)"
	echo "  -b database to search in options are GenBank or RefSeq"
	echo "  -h  display this help"
	echo "  -d  target download directory"
	echo
	echo "Author: Jake Hudson"
	echo
}

#################################
# Setup
#################################
links_file="tmp_ftp_links.txt"


#################################
# Process input options
#################################

# Get Options
while getopts ":hs:d:b:" option;
do
	case $option in
		h) # Display help
		  Help
		  exit;;
		s) # Search term for NCBI
		  search=$OPTARG;;
		d) # Directory to download into
		  dir=$OPTARG;;
		b) #Search genbank or refset
		  db=$OPTARG;;
		\?) #Invalid option
		  echo "Error: Invalid option"
		  Help
		  exit;;
		esac
done

###############################
# Main Program
###############################
# Make download directory 
echo "Downloading into $dir"
mkdir -p $dir

#get current date for log file
today=$(date -I | sed 's/-//g')
echo "Today is $today"

#search using NCBI eutils
echo ""
echo "searching for "${search}" in FtpPath_$db"
echo "search term is"
echo 'esearch -db assembly -query "${search}" |
		efetch -format docsum |
 			xtract -pattern DocumentSummary -element FtpPath_$db |
				sed 's/^ftp//' > $links_file'

esearch -db assembly -query "${search}" |
		efetch -format docsum |
 			xtract -pattern DocumentSummary -element FtpPath_$db |
				sed 's/^ftp//' > $links_file

# Give me esearch output 
echo ""
echo "esearch says:"
esearch -db assembly -query $search

# Count how many genome file I already have
ass_count=$(find $dir -type f -name *_genomic.fna.gz | wc -l)

# Count hoe may assembly reports I already have
rep_count=$(find $dir -type f -name *_assembly_report.txt | wc -l)

# Confirm if we want to proceede 
count=$(cat ${links_file} | wc -l)

echo ""
echo "You already have $ass_count assemblies downloaded"
echo "You already have $rep_count reports downloaded"
echo "" 

while true
	do
		read -p "${count} FtpPath_${db} links found, continue? y/n: " yn
		case $yn in 
			[Yy]* ) # Proceede with download
			          echo "Downloading now"
			          break;;
			[Nn]* ) # Abort Download
			          echo "Aborting"
			          rm ${links_file}
			          rmdir ${dir}
			          exit;;
			* ) #Invalid options
			    echo "Please answer y or n"
		esac
done

#Load list of ftp links
Links=$(cat $links_file)

# See if we have already completed files
comp_dls=$(ls ${dir})

#make header for rsync log file
logfile=${today}_sequence_downloader_log_rep.tsv
echo -e "genome\tftp_link\ttime_downloaded\tstatus" > $logfile

#make a counter so we can count progress done
counter=1

echo "starting downloading"
for ftp_address in $Links
do
	# Add some bits to the link so it only downloads the genomic fasta file and report
	accession=${ftp_address##*/}
	genome_link=${ftp_address}/${accession}_genomic.fna.gz
	report_link=${ftp_address}/${accession}_assembly_report.txt
	
	# Make directory for the files
	mkdir -p ${dir}/${accession}

#	echo ""
#	echo "ftp_address = ${ftp_address}"
#	echo "accession   = ${accession}"
#	echo "link        = ${Link}"
#	echo ""

	#get current time
	t=$(date +"%T")	

	#check if file has already been downloaded so we don't need to try again
#	if [[ $comp_dls =~ (^|[[:space:]])"$genome"($|[[:space:]]) ]]
#	then
#		
#		echo -e "${acession}\t${Link}\t${t}\tallready_downloaded" >> $logfile
#
#	else
#		#run rsync, but also write what has failed 
#		rsync --copy-links --times --recursive "rsync${genome_link}" ${dir} && \
#			echo -e "${genome}\t${Link}\t${t}\tsuccess" >> $logfile || \
#			echo -e "${genome}\t${Link}\t${t}\tfail" >> $logfile
#		
#		#0.4s sleep as NCBI has a limit of >3 access requests per second
#		sleep 0.4s
#	fi

	# Check if we do not have the genome file already
	if [[ ! -f ${dir}/${accession}/${accession}_genomic.fna.gz ]] 
	  then 
	    echo "downloading ${accession}_genomic.fna.gz"
	    rsync --times "rsync${genome_link}" ${dir}/${accession}/
	    sleep 0.4s
	fi

	# Check if we do not have the genome report already
	if [[ ! -f ${dir}/${accession}/${accession}_assembly_report.txt ]] 
	  then 
            echo "downloading ${accession}_assembly_report.txt"
	    rsync --times "rsync${report_link}" ${dir}/${accession}/
	    sleep 0.4s
	fi 
	  

	#calculate percentage complete
	per_comp=$(echo "scale=2 ; $counter * 100 / $count" | bc)

	echo "${per_comp}% complete"

	#add 1 to counter
	let counter=counter+1

done

