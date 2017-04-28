#!/bin/bash
# This script counts the number of lines in several files and put the values 
# into a table. It takes 2 arguments: a list of files and a type of gene
# sets. Here, it is used to make table summarizing the number of genes in
# each host group.

# Reinitializing table
echo '' > $2

# Looping over list of files
for i in $1;
do
    j=basename $i # extracting filename without path
    echo ${j%_*}','$(cat $i | wc -l) >> $2; 
    # adding a line to the table in the form "host_group","number of genes"
done
