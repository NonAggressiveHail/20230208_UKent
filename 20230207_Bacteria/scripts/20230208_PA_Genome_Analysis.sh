#script to analyse the genomes of P. aeruginosa for unknown gene content
#goal is to produce a graph which shows % of unknown genes
#uses conda environment "orf_prediction"

#setup what I am running on 
org="PA01"
dir="../data/geneomes/"

#find geneome file
genome=$(find ${dir}${org} -type f -name "*_genomic.fna")
echo ${genome}

#get accession number
tmp=${genome%_genomic.*}
accession=${tmp##*/}

#make protein output filename
gene_out=${dir}${org}/${accession}_genes.fna

#make genome output filename
prot_out=${dir}${org}/${accession}_proteins.faa

#run prodigal to predict genes and proteins
prodigal -i $genome -o $gene_out -a $prot_out
