
# The purpose of this script is to extract a list of genes that are orthologous  between bumblebee (B)- and honeybee (H)-specific
# bacterial species. It takes an ortholog table as input and output a similar table containing only B-H orthologs.
# Laurent Casini, Cyril Matthey-Doret
# 04.04.2017

#==========================================
## Loading data

# List of genomes IDs
ID_H <- c() # Honeybee-bacterial genomes
ID_B <- c() # Honeybee-bacterial genomes
ID_O <- c() # Honeybee-bacterial genomes

# Loading ortholog table line by line
ortho_tab = open("Genefamilies_all.txt",'r')
#==========================================

for line in ortho_tab:
    
