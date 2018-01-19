# This script extracts the GO terms from the eggNOG annotation output, counts
# the number of occurences of each term and performs a GO term enrichment test
# for a given host group (host group input file given as command line argument).
# The enrichment test uses the Fisher exact test and p-values are corrected for
# multiple testing using the Benjamini-Hochberg procedure. The script finally
# writes a list of significantly enriched (q-value <= 0.05 ) GO terms to a file,
# with their GO informations, oddratios, p-values and q-values.
# Cyril Matthey-Doret
# 18.05.2017

import re  # Regular expressions support
from sys import argv  # Using command line arguments
import pandas as pd  # Convenient DataFrames
#import MySQLdb  # Sending SQL queries to GO consortium database
import mysql.connector
from os.path import basename  # Allows to remove path from filename
from numpy import median  # fast vectorized median method
from copy import copy  # copy data structures
from scipy.stats import fisher_exact  # Used for enrichment test
import numpy as np  # Fast data structures and vectorized methods

""" Loading data """

# group_file = open('../data/gene_sets/BH_genes.txt','r')
group_file = open(argv[1],'r')  # OrthoMCL table for input host group
host_group = basename(argv[1]).split('_')[0]  # Extracting host group name
group_content = group_file.read()  # Reading whole file as a string

# Splitting strain|gene doublets
group_content = re.split(r'[\n\t]',group_content)
group_content = filter(None, group_content)
# In case there are consecutive tabs/newlines, empty items are removed

if len(argv) == 3:  # If a second file was given, compare first vs second
    group2_file = open(argv[2],'r')
    group2_content = group2_file.read()
    host_group2 = basename(argv[2]).split('_')[0]
    group2_content = re.split(r'[\n\t]',group2_content)
    group2_content = filter(None, group2_content)

eggNOG = pd.read_csv('data/eggNOG_all_annotations.txt',sep='\t',header=None)
# Reading eggNOG annotations (concatenated from all other strains)

""" Processing annotations """

GO = eggNOG[[0,5]]  # Only keeping GO terms and corresponding gene ID

GO_group = GO[GO[0].isin(group_content)].dropna()
# Keeping only Genes that are in orthoMCL table of current host group
if len(argv) == 3:  # For dual group comparison
    GO_set = GO[GO[0].isin(group2_content)].dropna()
else:  # For one versus all others comparison
    GO_set = copy(GO).dropna()  # For all strains Duplicating table and removing NAs
    # Excluding rows in host group to keep only all other strains (not in host group)
    GO_set = GO_set.drop(GO_group.index)

# Splitting GO terms into a list for each gene
GO_group[5] = GO_group[5].str.split(',')  # GO terms in host group
GO_set[5] = GO_set[5].str.split(',')  # GO terms in all other strains

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
        go_dict[go] = 1  # If absent, add it as a new key

GO_occ, GO_all = {},{}  # Inializing dict. for host group and all other strains
GO_group[5].map(lambda x: map(lambda y: GO_expand(y,GO_occ),x))  # host group
GO_set[5].map(lambda x: map(lambda y: GO_expand(y,GO_all),x))  # other strains
# Using nested maps to apply 'GO_expand' on each GO term.
# The first map returns each cell of column '5' and the second map will apply
# the function of every element of the list in each cell.

"""
Fisher exact test
table:
    | host | all |
 GO | GOh  | GOa |
!GO | nGOh | nGOa|

Question:
    Is a GO term more frequent in the 'host' group than in 'all other' strains?
general concept:
    if Goh/nGOh >> GOa/nGOa -> GO enriched in host group
    This is given by the "odd ratio": OR = (GOh*nGOa)/(GOa*nGOh)
implementation:
    fisher_exact([[GOh,GOa],[nGOh,nGOa]])
"""

# Building dataframe with all GO terms found in host group.
# GOh: number of occurences of each GO term in host group
# GOa: number of occurences of corresponding GO terms in all other strains

GO_calc = pd.Series(GO_occ).to_frame(name='GOh').merge(
    pd.Series(GO_all).to_frame('GOa'),how='left',right_index=True,
    left_index=True)

GO_calc.GOa = GO_calc.GOa.fillna(value=0)  # Replacing NA's with 0's

# Storing total number of GO terms occuring in host group and all other strains
tot_GO = {'host':sum(GO_occ.values()),
          'all':sum(GO_all.values())}

# adding column for all occurences of other GO terms in:
GO_calc['nGOh'] = GO_calc.GOh.apply(lambda x: tot_GO['host'] - x)  # host group
GO_calc['nGOa'] = GO_calc.GOa.apply(lambda x: tot_GO['all'] - x)  # all other

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


ebi_connect = mysql.connector.connect(host = "mysql.ebi.ac.uk",user="go_select",
                             passwd="amigo",database="go_latest",port=4085)
# Establishing connection with EBI mirror of GO database

query = "SELECT * FROM term WHERE acc IN ('%s')" % "','".join(map(str,list(
    GO_enrich.index)))
# SQL query string requesting all db entries corresponding to current host group
print('Sending SQL request for GO terms of group ' + host_group +'...' )
sql_annot = pd.read_sql(query, con=ebi_connect)
# Saving server response as dataframe

ebi_connect.close()  # Closing connection to SQL server

""" Filtering annotations """

full_annot = sql_annot.merge(GO_enrich, left_on='acc',
                             right_index=True,how='inner')
# Merging frequency and annotations into single dataframe

top_GO = full_annot.loc[(full_annot['qval']<=np.float(0.05))]
# Including only annotations with q-values below cutoff
if len(argv) == 3:  # Writing dataframe to output file (group1 vs group2)
    top_GO.to_csv('data/annotations/'+host_group+'vs'+host_group2+'_annot.txt',
                  sep='\t')
else:  # Writing dataframe to output file (group versus all others)
    top_GO.to_csv('data/annotations/'+host_group+'_annot.txt',sep='\t')

print('Enrichment analysis of group ' + host_group +' performed successfully !' )
