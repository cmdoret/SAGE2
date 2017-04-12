#!/bin/bash

cd ./data/gene_sets/
echo '' > gene_number.csv

for i in *.txt;
do
    echo ${i%_*}','$(cat $i | wc -l) >> gene_number.csv;
done
