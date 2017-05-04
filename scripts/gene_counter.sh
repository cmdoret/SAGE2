#!/bin/bash
# This script counts the number of lines in several files and put the values 
# into a table. It takes 2 arguments: the table filename and a common pattern
# in the filename of all sets of interest. Here, it is used to make table 
# summarizing the number of genes in each host group.

# Reinitializing table
echo -n '' > $1

# Looping over list of files
for i in data/gene_sets/*$2;
do
    j=$(basename $i) # extracting filename without path
    echo ${j%_*}','$(cat $i | wc -l) >> $1; 
    # adding a line to the table in the form "host_group","number of genes"
done
