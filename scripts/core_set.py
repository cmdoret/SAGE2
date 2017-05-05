from sys import argv
from os.path import join
from copy import copy
import re
# This script extracts gene families common to all strains in host groups given
# as command line arguments B (Bumblebee), H (Honeybee) and O (Outgroup).
# Laurent Casini, Cyril Matthey-Doret
# 25.04.2017

# Only keeping parameters containing one species group
param = [h for h in argv[1:] if re.compile(r'[BHO]').search(h)]
groups = []

# Parsing parameters so that BHO works like B H O
for p in param:
    for c in p:
        groups += c

# Preventing duplicate groups
group_set = set(groups)

# List of all strains:
all_strains = {'H':['JF72','F259','JG30','JF76','F260','F261','F262','F263','JF73','JF74',
     'JF75','WB8','WB10','L185','L186','L184','L183'],
'B':['F225','F230','F233','F234','F236','F237','F228','F245','F246','F247'],
'O':['LA14','LA2','LDB','LGAS','LHV','LJP','WANG','JG29']}

# Subsetting strains with hosts matching input parameters
strains=[]
for g in group_set:
    strains += all_strains[g]

# Generating filenames from input parameters
in_name = ''.join(sorted(group_set)) + "_genes.txt"
out_name = ''.join(sorted(group_set)) + "_core_set.txt"

outfile=open(join("data","gene_sets",out_name),'w') # Creating new file/erasing previous version
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
