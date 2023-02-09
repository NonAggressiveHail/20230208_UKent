#!bin/bash
#script to download all complete and chromosome level pseudomonas refseq assemblies from NCBI
#author: Jacob Hudson
#NB can be run on login node
#NB reccomend running without downloading anything first, to check how many files will be downloaded
#activate conda
source activate eutils

#set wd
cd ../

#get current date for log file
today=$(date -I | sed 's/-//g')

#search using NCBI eutils
esearch \
	-db assembly \
	-query "Pseudomonas [ORGN] AND latest [filter]" |
		efetch -format docsum |
			xtract -pattern DocumentSummary -element FtpPath_GenBank |
				sed 's/^ftp//' > ./paper_data/all_pseudo_spp_genbank_ftp_links.txt


#print how many genomes will be downloaded 
count=$(cat ./paper_data/all_pseudo_spp_genbank_ftp_links.txt | wc -l)
echo "downloading ${count} genomes"

##Load list of ftp links
#Links_File="./paper_data/all_pseudo_spp_refseq_ftp_links.txt"
#Links=$(cat $Links_File)
#
##load completed files
#comp_dls=$(ls ./raw_data/all_pseudo_spp_refseq/)
#
##make header for rsync log file
#logfile=./scripts/${today}_sequence_downloader_log_rep_${i}.tsv
#echo -e "genome\tftp_link\ttime_downloaded\tstatus" > $logfile
#
##make a counter so we can count progress done
#counter=1
#
#for Link in $Links
#do
#
#	#get name of genome I am downloading
#	genome=${Link##*/}
#	
#	#get current time
#	t=$(date +"%T")	
#
#	#check if file has already been downloaded so we don't need to try again
#	if [[ $comp_dls =~ (^|[[:space:]])"$genome"($|[[:space:]]) ]]
#	then
#		
#		echo -e "${genome}\t${Link}\t${t}\tallready_downloaded" >> $logfile
#
#	else
#		#run rsync, but also write what has failed 
#		rsync --copy-links --times --verbose --recursive "rsync${Link}" ./raw_data/all_pseudo_spp_refseq/ && \
#			echo -e "${genome}\t${Link}\t${t}\tsuccess" >> $logfile || \
#			echo -e "${genome}\t${Link}\t${t}\tfail" >> $logfile
#		
#		#0.4s sleep as NCBI has a limit of >3 access requests per second
#		sleep 0.4s
#	fi
#
#
#	#calculate percentage complete
#	per_comp=$(echo "scale=2 ; $counter * 100 / $count" | bc)
#
#	echo "${per_comp}% complete"
#
#	#add 1 to counter
#	let counter=counter+1
#
#done
