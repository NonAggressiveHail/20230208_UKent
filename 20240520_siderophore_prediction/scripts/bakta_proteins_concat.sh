#!/bin/bash
# Script which combines together the "proteins" and "hyptothetical proteins" output fasta files from bakta
# Needed for orthofinder to use 

for dir in ../data/bakta/*
do
	
	name=$(basename $dir)
		
	cat ${dir}/${name}.faa ${dir}/${name}.hypotheticals.faa > ${dir}/${name}_all_genes.faa
	#echo "Command: cat ${dir}/${name}.faa ${dir}/${name}.hypotheticals.faa > ${dir}/${name}_all_genes.faa"
	echo "Done ${name}"
	
done