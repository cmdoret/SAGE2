from time import sleep
from copy import copy
# This script extracts genes common to all strains.
# Laurent Casini, Cyril Matthey-Doret
# 25.04.2017

# List of strains:


strains = ['JF72','F259','JG30','JF76','F260','F261','F262','F263',
        'JF73','JF74','JF75','WB8','WB10','L185','L186',
        'L184','L183','F225','F230','F233','F234','F236','F237',
        'F228','F245','F246','F247','LA14','LA2','LDB','LGAS','LHV','LJP','WANG',
        'JG29']

outfile=open("core_set.txt",'w')
k=0
with open("data/gene_sets/BHO_genes.txt") as ortho:
    for line in ortho:
        tmp_str = copy(strains)
        gene = line.split("\t")
        for g in gene:
            for s in tmp_str:
                if g.split("|")[0] == s:
                    tmp_str.remove(s)
        if len(tmp_str) == 0:
            outfile.write(line)
        else:
            #print(len([x for x in strains if x not in tmp_str]))
            k+=1

outfile.close()
print(k)
