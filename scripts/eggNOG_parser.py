# This script extracts the GO terms from the eggNOG annotation output and counts
# produces a list of GO terms with their frequency in a given host group.
# Cyril Matthey-Doret
# 18.05.2017
# TODO: REPLACE PATHS TO MAKEFILE COMPATIBLE VERSIONS

import re  # Regular expressions support
from sys import argv  # Using command line arguments
import pandas as pd  # Convenient DataFrames
import MySQLdb  # Sending SQL queries to GO consortium database
from os.path import basename
from numpy import median
from copy import copy



""" Loading data """

#group_file = open('data/gene_sets/BH_genes.txt','r')
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

# Building GO occurences dictionary from dataframe
GO_occ = pd.Series(GO_occ)  # Restructuring dict as Series
GO_calc = GO_occ.to_frame(name='occ').merge(pd.Series(GO_all).to_frame('all'),
                                              how='left',right_index=True,
                                              left_index=True)

GO_freq = GO_calc['occ']/GO_calc['all']
# Frequencies are computed as occurences in host group/occurences in all strains

""" Querying database """

ebi_connect = MySQLdb.connect(host = "mysql.ebi.ac.uk",user="go_select",
                             passwd="amigo",db="go_latest",port=4085)
# Establishing connection with EBI mirror of GO database

query = "SELECT * FROM term WHERE acc IN ('%s')" % "','".join(map(str,list(
    GO_freq.index)))
# SQL query string requesting all db entries corresponding to current host group

sql_annot = pd.read_sql(query, con=ebi_connect)
# Saving server response as dataframe

ebi_connect.close()  # Closing connection

""" Filtering annotations """

full_annot = sql_annot.merge(GO_freq.to_frame(name='freq'), left_on='acc',
                             right_index=True,how='inner')
# Merging frequency and annotations into single dataframe

full_annot = full_annot[full_annot.freq > full_annot.freq.median()]
# Including only annotations with more than median frequency
annot_bytype = full_annot.groupby('term_type')
# Group annotations by category (mol. function, bio. process or cell. component)
top_GO=annot_bytype.apply(lambda x: x.sort_values('freq',ascending=False)[0:5])
# Sort annotations by number of occurences and only selecting 5 top annotations
# for each GO term category

host_group = basename(argv[1]).split('_')[0]
top_GO.to_csv('data/annotations/'+host_group+'_annot.txt',sep='\t')
# Writing dataframe to output file
