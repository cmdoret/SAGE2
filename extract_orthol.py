
# The purpose of this script is to extract a list of genes that are orthologous  between bumblebee (B)- and honeybee (H)-specific
# bacterial species. It takes an ortholog table as input and output a similar table containing only B-H orthologs.
# Laurent Casini, Cyril Matthey-Doret
# 04.04.2017

#==========================================
## Loading data

# List of genomes IDs
ID_H = ["F262","F263","F261"] # Honeybee-bacterial genomes
ID_B = ["F225","F228","F233","F237"] # Bumblebee-bacterial genomes
ID_O = [] # Outgroup bacterial genomes

# Loading ortholog table and creating output files
ortho_tab = open("../Genefamilies_all.txt",'r')
bumble_genes = open("bumble_genes.txt",'w')
honey_genes = open("honey_genes.txt",'w')
#==========================================

for line in ortho_tab:  # Each line is a gene family
    tmp = line.split("\t")  # split line by word (tab separated)
    B = False; H = False  # Is this family present in bumble/honeybees
    for word in tmp: # For each word (a word is: "strain|gene")
        gene = word.split("|")  # sparates strain from gene in a list
        if gene[0] in ID_H: # if the strain is in honeybee
            H = True
        if gene[0] in ID_B: # if the strain is in bumblebee
            B = True
    if B and not H:  # if gene family only present in bumblebees
        bumble_genes.write(line)  # writing family to bumblebee file
    if H and not B:  # if gene family only present in honeybees
        honey_genes.write(line)  # writing family to honeybee file
