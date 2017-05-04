import itertools
from os.path import join
import pandas as pd
# The purpose of this script is to extract a list of genes that are orthologous
# between bumblebee (B)- and honeybee (H)-specific bacterial species. It takes
# an ortholog table as input and output a similar table containing only B-H
# orthologs.
# Laurent Casini, Cyril Matthey-Doret
# 04.04.2017

#==========================================
## Loading data
files = {}  # Data structure containing output files
groups = ['B','H','O']
# Three groups considered here: Bumblebee , Honeybee, Outgroup
for L in range(1,len(groups)+1):
    for combi in itertools.combinations(groups,L):
        # Need all combinations of 1, 2 or all 3 of these groups
        combi = ''.join(combi)  # Tuple output concatenated to string
        files[combi] = open(join("./data","gene_sets",combi + "_genes.txt"),'w')
        # Each file is a dictionary value with the group as key (e.g. "BH").
gnm = pd.read_csv(join("./data","strain_list"),sep="\t")
# Contains genomes ID's
ID = {}
for g in groups:
    ID[g] = list(gnm[gnm.host_group==g].OrthoMCL_prefix)
# Honeybee, bumblebee and outgroup-bacterial genomes

# Loading ortholog table and creating output files
ortho_tab = open("./data/mclOutput",'r')
#==========================================

for line in ortho_tab:  # Each line is a gene family
    tmp = line.split("\t")  # split line by word (tab separated)
    comb_group = set()  # Is this family present in bumble/honeybees/outgroup
    for word in tmp: # For each word (a word is: "strain|gene")
        gene = word.split("|")  # sparates strain from gene in a list
        for g in groups:
            if gene[0] in ID[g]: # if the strain is in ...
                comb_group.add(g)
    comb_group = ''.join(sorted(comb_group))  # compressing set to a string
    files[comb_group].write(line)  # Writing line to file matching group
ortho_tab.close()  # closing input file
for x in files: files[x].close()  # closing all output files in a loop
