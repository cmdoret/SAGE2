
# To switch between the complete and reduced sets of strains, uncomment the desired
# GFAM and SLIST files and comment the other

# Complete gene family table (GFAM) and strain list (SLIST)
GFAM=./data/mclOutput.txt
SLIST=./data/strain_list.tsv

# Reduced gene family table (GFAM) and strain list (SLIST)
#GFAM=./data/mclOutput_GFs_reduced.txt
#SLIST=./data/strain_list_reduced.tsv

# Path to output folders
SDIR=./data/gene_sets/
ANNOT=./data/annotations
FREQ=./data/frequencies

# Host groups to consider (i.e. all combinations of B, H and O by default)
GROUPS=B BH BHO BO H HO O
