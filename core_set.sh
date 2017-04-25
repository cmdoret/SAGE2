#!/bin/bash 

# This script extracts genes common to all strains.
# Laurent Casini, Cyril Matthey-Doret
# 25.04.2017

# List of strains:

strains=(JF72 F259 JG30 JF76 F260 F261 F262 F263 JF73 JF74 JF75 WB8 WB10 L185 L186 L184 L183 F225 F230 F233 F234 F236 F237 F228 F245 F246 F247 LA14 LA2L DBLG ASLH VLJP WANG JG29)

echo -n "" > core_set.txt

while read p; do

    strain_line=$(echo $p | awk 'BEGIN {RS="\t";FS="|";x=""}
    {print $1}')
    #echo $strain_line
    sleep 1
done <data/gene_sets/BHO_genes.txt
