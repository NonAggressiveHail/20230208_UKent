#!/bin/bash
# Count number of files in a folder

for folder in ../raw_data/P_aeruginosa_genomes/*
do

 n_contents=$(ls $folder | wc -l)

 if [[ $n_contents != 2 ]]
   then
     echo "${folder} has ${n_contents}"

 fi
done


