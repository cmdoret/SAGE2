# This script extracts the GO terms from the eggNOG annotation output and counts
# produces a list of GO terms with their frequency in a given host group.
# Cyril Matthey-Doret
# 18.05.2017
# TODO: REPLACE PATHS TO MAKEFILE COMPATIBLE VERSIONS

import re  # Regular expressions support
from sys import argv  # Using command line arguments
import pandas as pd  # Convenient DataFrames
import MySQLdb  # Sending SQL queries to GO consortium database



""" Loading data """

group_file = open('../data/gene_sets/BH_genes.txt','r')
#group_file = open(argv[1],'r')  # OrthoMCL table for input host group
group_content = group_file.read()  # Reading whole file as a string
group_content = re.split(r'[\n\t]',group_content)
# Splitting strain|gene doublets
eggNOG = pd.read_csv('../data/eggNOG_all_annotations',sep='\t',header=None)
# Reading eggNOG annotations (concatenated from all strains)
strains = pd.read_csv('../data/strain_list',sep='\t')
# Strain list for ID conversions

""" Processing annotations """

GO = eggNOG[[0,5]]  # Only keeping GO terms and corresponding gene ID
# GO[[0]].applymap(lambda x: x.replace('_','|'))

group_GO = GO[GO[0].isin(group_content)].dropna()
# Keeping only Genes that are in orthoMCL table of current host group
group_GO[5] = group_GO[5].str.split(',')
# Splitting GO terms into a list for each gene


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

GO_freq = {}
group_GO[5].map(lambda x: map(lambda y: GO_expand(y,GO_freq),x))
# Building GO occurences dictionary from dataframe
GO_freq = pd.Series(GO_freq)
top_GO = GO_freq.sort_values(ascending=False)[:10]

""" Querying database """

connection = MySQLdb.connect(host = "mysql.ebi.ac.uk",user="go_select",
                             passwd="amigo",db="go_latest",port=4085)
# Establishing connection with EBI mirror of GO database

cursor = connection.cursor()  # Emulating cursor
cursor.execute ("SELECT * FROM term WHERE acc='GO:0005634'")
# Sending query for current term
row = cursor.fetchone()  # Return first row of server response
cursor.close()  # Removing cursor
connection.close()  # Closing connection
