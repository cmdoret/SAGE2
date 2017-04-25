from copy import copy
# This script extracts gene families common to all strains.
# Laurent Casini, Cyril Matthey-Doret
# 25.04.2017

# List of all strains:
strains = ['JF72','F259','JG30','JF76','F260','F261','F262','F263',
        'JF73','JF74','JF75','WB8','WB10','L185','L186',
        'L184','L183','F225','F230','F233','F234','F236','F237',
        'F228','F245','F246','F247','LA14','LA2','LDB','LGAS','LHV','LJP','WANG',
        'JG29']

outfile=open("data/core_set.txt",'w') # Creating new file/erasing previous version
with open("data/gene_sets/BHO_genes.txt") as ortho:
    # Opening Mcl ortholog table
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
