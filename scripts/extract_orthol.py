
# The purpose of this script is to extract a list of genes that are orthologous
# between bumblebee (B)- and honeybee (H)-specific bacterial species. It takes
# an ortholog table as input and output a similar table containing only B-H
# orthologs. The second part of the script appends genes uniques to a single
# strain to the output table. It uses the fasta files of each strain to extract
# the genes.
# Laurent Casini, Cyril Matthey-Doret
# 04.04.2017

import itertools
from os.path import join
import pandas as pd
from Bio import SeqIO
import re

## Loading data
files = {}  # Data structure containing output files
groups = ['B','H','O']
# Three groups considered here: Bumblebee , Honeybee, Outgroup
for L in range(1,len(groups)+1):
    for combi in itertools.combinations(groups,L):
        # Need all combinations of 1, 2 or all 3 of these groups
        combi = ''.join(combi)  # Tuple output concatenated to string
        files[combi] = open(join("data","gene_sets",combi + "_genes.txt"),'w')
        # Each file is a dictionary value with the group as key (e.g. "BH").
gnm = pd.read_csv(join("../data","strain_list"),sep="\t")
# Contains genomes ID's
ID = {}
for g in groups:
    ID[g] = list(gnm[gnm.host_group==g].OrthoMCL_prefix)
# Honeybee, bumblebee and outgroup-bacterial genomes

# Loading ortholog table and creating output files
ortho_tab = open(join("data","mclOutput"),'r')

## Parsing MCL table
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

## Append unique genes to each table
fasta = {} # List of each strain with its parsed FASTA content
# Structure of 'fasta': {'H':{'strain1':['gene1','gene2']},...}
for group in ID:  # Iterating over hosts (H, B and O)
    fasta[group] = []  # Initiating new dictionary for FASTA parsing
    for strain in ID[group]:  # Iterating over strains in given group
        prefix = gnm[gnm.OrthoMCL_prefix==strain].Strain_prefix.values[0]
        # Extracting strain prefix from strain list for given orthoMCL prefix
        fa_file = join('data','genomes',prefix + '.faa')
        fasta[group] += [strain + '|' + filter(None,re.split("[|_]",gene.id))[-1]
                               for gene in SeqIO.parse(fa_file,'fasta')]
        # Using BioPython to parse fasta file and store all genes in dictionary
        # Only the gene ID is extracted and the strain name is added to form
        # an orthoMCL-compatible entry.

    uniq_list = []  # List to store unique genes.
    gr_ortho = open(join('data','mclOutput'),'r').read()
    #gr_ortho = open(join('data','gene_sets',group + '_genes.txt'),'r').read()
    # Storing content of orthoMCL file for current group
    track_uniq = [0,0]  # Tracking [total, unique] genes for info
    for gene in fasta[group]:
        # Iterating over genes from fasta entries
        if gene in gr_ortho:  # Checking if gene is already in ortho table
            track_uniq[0] += 1  # incrementing total gene count
        else:
            track_uniq[1] += 1  # Incrementing unique gene count
            uniq_list.append(gene)  # If not, unique gene -> append to list
    print(('{0} genes added to orthoMCL table of group {1} out of {2} in FASTA'+
           ' files.').format(track_uniq[1],group,track_uniq[0]+track_uniq[1]))
    add_ortho = open(join('data','gene_sets',group + '_genes.txt'),'a')
    # Opening file in append mode to write unique genes at the end
    for uniq in uniq_list: add_ortho.write(uniq +'\n')
