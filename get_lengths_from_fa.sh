#!/bin/bash

#Change seq-dir to be the location of the fasta's used to get abundance (probably from CoPTR unless we do something different)
#stats.tsv is created in the directory from which this is run, use it in the python script
#ASSUMES THAT ALL FILES IN DIRECTORY ARE FASTA

#create a tsv using seqkit stats
#note, -j indicates CPU number
seqkit stats --infile-list <(find seq-dir/ -name "*.*") --tabular -j 12 -o stats.tsv

#remove the file extensions with an awk substitution
awk -F'\t' '{sub(/\.[^.]+$/, "", $1)} 1' stats.tsv > temp.tsv && mv temp.tsv stats.tsv

#the first line creates space delimited instead of tab delimeted, replace spaces with tabs
sed 's/ /\t/g' stats.tsv > temp.tsv && mv temp.tsv stats.tsv
