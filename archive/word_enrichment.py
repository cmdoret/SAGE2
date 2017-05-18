# This is an exploratory script! It extracts most frequent words in the whole
# annotation table and us it to compute word enrichment in an orthology table.
# Cyril Matthey-Doret
# 12.05.2017

# TODO: Find workaround to initialize dictionnary without using a set
# TODO: Use pandas methods to remove top_group entries that are in exclude list
# TODO: Compute ratio of group/overall freq for each top word to get enrichment

from sys import argv
import pandas as pd
import re
import pickle

loc_interest = ["nucleus","membrane","cytoplasm","extracell","mitochondri"]
func_interest = ["replication","repair","transduction","transport","metal",
                 "sugar degradation","tRNA","ribosom","ATP","histone"]
exclude = ["peg","(EC","hypothetical","protein","component","ff","","cluster",
           "of","with","in","specific","family","bacterial","binding","and",
           "(TC","subunit"]

def word_counter(input_string):
    """
    This function counts the number of occurences of each word in a string.

    :param input_string: input string in which to count words.
    :returns: a dictionary of words and their respective frequencies
    """
    splitted = re.split(r'[_\-\n\t ]',input_string)
    # Splitting string by several different separators using regex
    word_set = set(splitted)  # Only keeping unique words
    word_freq = {word:0 for word in word_set}  # Initializing dictionary
    pw = ''; n = 0  # n: keeps track of current word occurences; pw: previous word
    for w in sorted(splitted):
        # Iterating on word (in alphabetically sorted list)
        if pw == w:  # If word did not change
            n += 1  # Add 1 occurence
        else:
            if pw:  # Prevents adding empty string to the dict (first iteration)
                word_freq[pw] += n  # Writing number of occurences of current word
            n = 1  # Resetting number of occurences (1 observation, caused change)
            pw = w  # Updating previous word
    return word_freq

try:
    pickle.load(open('all_words_freq','r+b'))
    # Trying to load pickled dictionary to avoid recalculating word frequencies
except FileNotFoundError:
    # If file does not exist, counting words and storing pickled dictionary
    whole_table = open('../data/all_strains_annot.tsv','r')
    whole_content = whole_table.read()
    word_freq = word_counter(whole_content)
    # Counting frequency of every word
    pickle.dump(word_freq,open('all_words_freq','w+b'))
    # Saving pickled dictionary for later runs
    whole_table.close()  # Closing file connection

annot = open(argv[1],'r')  # Input annotation file for given group
group_words = word_counter(annot.read())
top_group = pd.Series(group_words)
val_group = [f for f in top_group if f not in exclude]
print(top_group[val_group].sort_values(ascending=False)[:20])
