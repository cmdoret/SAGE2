
# This script extracts gene families common to all strains in host groups given
# as command line arguments B (Bumblebee), H (Honeybee) and O (Outgroup).
# Laurent Casini, Cyril Matthey-Doret
# 25.04.2017

from sys import argv
from os.path import join
from copy import copy
import pandas as pd
import re

s_path = argv[1]  # path to list of strains
# Only keeping parameters containing one species group
param = [h for h in argv[2:] if re.compile(r'[BHO]').search(h)]
groups = []

# Parsing parameters so that BHO works like B H O
for p in param:
    for c in p:
        groups += c

# Preventing duplicate groups
group_set = set(groups)

# List of all strains:
s_tbl = pd.read_csv(s_path, sep='\t')
all_strains = {k:list(s_tbl.loc[s_tbl.host_group == k,"OrthoMCL_prefix"])
               for k in ['B','H','O']}


# Subsetting strains with hosts matching input parameters
strains=[]
for g in group_set:
    strains += all_strains[g]

# Generating filenames from input parameters
in_name = ''.join(sorted(group_set)) + "_genes.txt"
out_name = ''.join(sorted(group_set)) + "_core_set.txt"

outfile=open(join("data","gene_sets",out_name),'w')
# Creating new file/erasing previous version
with open(join("data","gene_sets",in_name),'r') as ortho:
    # Opening MCL ortholog table
    for line in ortho:
        # iterating over lines of ortholog table (gene families)
        tmp_str = copy(strains)
        # copying a temporary list of strains for each new line
        gene = line.split("\t")
        # splitting gene family into list of strain|gene pairs
        for g in gene:
            # iterating over pairs
            for s in tmp_str:
                # iterating over strains in temporary list
                if g.split("|")[0] == s:
                    # if strain of pair is found in list
                    tmp_str.remove(s)
                    # it is removed
        if len(tmp_str) == 0:
            # at the end of the line, checking if all strains were removed
            outfile.write(line)
            # if so, writing line to output file
outfile.close()
print("Extracted core gene families for group "+''.join(sorted(group_set)))
