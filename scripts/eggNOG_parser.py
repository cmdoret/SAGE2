# This script extracts the GO terms from the eggNOG annotation output and counts
# produces a list of GO terms with their frequency in a given host group.
# Cyril Matthey-Doret
# 18.05.2017

import re  # Regular expressions support
from sys import argv  # Using command line arguments
import pandas as pd  # Convenient DataFrames
import MySQLdb  # Sending SQL queries to GO consortium database
from os.path import basename  # Allows to remove path from filename
from numpy import median  # fast vectorized median method
from copy import copy  # copy data structures
from scipy.stats import fisher_exact
import numpy as np


""" Loading data """

# group_file = open('../data/gene_sets/BH_genes.txt','r')
group_file = open(argv[1],'r')  # OrthoMCL table for input host group
group_content = group_file.read()  # Reading whole file as a string
group_content = re.split(r'[\n\t]',group_content)
group_content = filter(None, group_content)
# Splitting strain|gene doublets

eggNOG = pd.read_csv('data/eggNOG_all_annotations',sep='\t',header=None)
# Reading eggNOG annotations (concatenated from all strains)
strains = pd.read_csv('data/strain_list',sep='\t')
# Strain list for ID conversions

""" Processing annotations """

GO = eggNOG[[0,5]]  # Only keeping GO terms and corresponding gene ID

GO_set = copy(GO).dropna()  # Duplicating table and removing NAs
GO_group = GO[GO[0].isin(group_content)].dropna()
# Keeping only Genes that are in orthoMCL table of current host group

# Splitting GO terms into a list for each gene
GO_group[5] = GO_group[5].str.split(',')  # GO terms in host group
GO_set[5] = GO_set[5].str.split(',')  # GO terms in all strains

def GO_expand(go, go_dict):
    """
    This function fills a dictionary with the number of occurences of GO terms.
    :param go: a GO term
    :param go_dict: a dictionary that will contain GO terms as keys and number
    of occurences as values
    :returns: None
    """
    try:
        go_dict[go] += 1  # If term already in dictionary, increment occurences
    except KeyError:
        go_dict[go] = 0  # If absent, add it as a new key

GO_occ, GO_all = {},{}  # Inializing dictionaries for host group and all strains
GO_group[5].map(lambda x: map(lambda y: GO_expand(y,GO_occ),x))  # host group
GO_set[5].map(lambda x: map(lambda y: GO_expand(y,GO_all),x))  # all strains
# Using nested maps to apply 'GO_expand' on each GO term.
# The first map will return each cell of column 5 and the second map will apply
# the function of every element of the list in each cell.

"""
Fisher exact test
table:
    | host | all |
 GO | GOh  | GOa |
!GO | nGOh | nGOa|

GO: Gene Ontology term
n: not
h: host
a: all

Question: Is a GO term more frequent in the 'host' group than in 'all' strains?
general concept: if Goh/nGOh >> GOa/nGOa -> GO enriched in group
implementation: fisher_exact([[GOh,GOa],[nGOh,nGOa]])
"""

# Building dataframe with all GO terms found in host group.
# GOh: number of occurences of each GO term in host group
# GOa: number of occurences of each GO term in all strains
GO_calc = pd.Series(GO_occ).to_frame(name='GOh').merge(
    pd.Series(GO_all).to_frame('GOa'),how='left',right_index=True,
    left_index=True)

# Storing total number of GO terms occuring in both host group and all strains
tot_GO = {'host':sum(GO_occ.itervalues()),
          'all':sum(GO_all.itervalues())}

# adding column for all occurences of other GO terms in:
GO_calc['nGOh'] = GO_calc.GOh.apply(lambda x: tot_GO['host'] - x)  # host group
GO_calc['nGOa'] = GO_calc.GOa.apply(lambda x: tot_GO['all'] - x)  # all strains

def Fisher_row(r):
    """
    This function takes a row of DataFrame as input and performs Fisher's exact
    test on it.
    :param r: pandas DataFrame row to be used as input for Fisher's exact test
    :return: a tuple containing the output of the test. Namely, the odds ratio
    and the p-value
    """
    row=[[r['GOh'],r['GOa']],[r['nGOh'],r['nGOa']]]
    return fisher_exact(row)

def p_adjust_bh(p):
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    :param p: a series of p-values
    :return: an array of BH-corrected p-values (q-values)
    """
    p = np.asfarray(p)  # Storing p-values as numpy array
    by_descend = p.argsort()[::-1]
    # getting array of indexes sorting p-values in descending order
    by_orig = by_descend.argsort()
    # getting array of indexes sorting by_descend indexes back (used to move
    # p-values in original order in the end)
    steps = float(len(p)) / np.arange(len(p), 0, -1)  # Storing increasing steps
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    # Correcting p-values and un-sorting back into original order
    return q[by_orig]

GO_calc['Fisher_out'] = GO_calc.apply(Fisher_row,axis=1)
# Use exact Fisher test to measure enrichment
GO_enrich=pd.DataFrame({'oddratios':GO_calc['Fisher_out'].apply(lambda r: r[0]),
                          'pval':GO_calc['Fisher_out'].apply(lambda r: r[1])})
# Splits Fisher test outputs into separate columns

GO_enrich['qval'] = p_adjust_bh(GO_enrich['pval'])  # Correct p-values with BH

""" Querying database """

ebi_connect = MySQLdb.connect(host = "mysql.ebi.ac.uk",user="go_select",
                             passwd="amigo",db="go_latest",port=4085)
# Establishing connection with EBI mirror of GO database

query = "SELECT * FROM term WHERE acc IN ('%s')" % "','".join(map(str,list(
    GO_enrich.index)))
# SQL query string requesting all db entries corresponding to current host group

sql_annot = pd.read_sql(query, con=ebi_connect)
# Saving server response as dataframe

ebi_connect.close()  # Closing connection to SQL server

""" Filtering annotations """

full_annot = sql_annot.merge(GO_enrich, left_on='acc',
                             right_index=True,how='inner')
# Merging frequency and annotations into single dataframe

top_GO = full_annot.loc[(full_annot['qval']<np.float(0.05))]
# Including only annotations with q-values below cutoff

host_group = basename(argv[1]).split('_')[0]
top_GO.to_csv('data/annotations/'+host_group+'_annot.txt',sep='\t')
# Writing dataframe to output file
