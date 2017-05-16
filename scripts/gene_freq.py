# This scripts counts the number of strains in which each gene family is present.
# Laurent Casini, Cyril Matthey-Doret
# 16.05.2017

from sys import argv  # Allows to use command line arguments
from os.path import basename,join  # Removes path from filename

ortho_file = argv[1]  # Path of input file
store_freq = []  # Initializing container for frequencies
with open(ortho_file,'r') as ortho:  # Opening file in read mode
    for line in ortho:  # Iterating over line
        genfam = line.split('\t')  # Splitting genes in each gene families
        store_freq.append(str(len(genfam)))
        # Appending number of occurences to store_freq

group_name = basename(ortho_file).split('_')[0]  # Extract group from filename
path = 'data/frequencies'
out_name = join(path,(group_name+'_gfreq.txt'))
out_freq = open(out_name,'w')  # Opening new output file for frequencies
print(','.join(store_freq))
out_freq.write(','.join(store_freq))  # Writing frequency list to output file
out_freq.close()
