import itertools
from os.path import join
import pandas as pd
# This script adds unique genes (i.e. genes present only in one strain) to the
# respective gene sets. It compares the FASTA file of each strain to the ortholog
# table of the corresponding group to those genes.

# Cyril Matthey-Doret
# 04.05.2017


for c in files:
    files[c].write('\n'.join(uniq_genes[c]))
