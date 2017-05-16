# This script produces a histogram given a series of numbers. In this context, 
# each number in the series correspond to a gene family and its value represent 
# the number of strains in which it is present.
# Cyril Matthey-Doret, Laurent Casini
# 16.05.2017

input_file <- commandArgs(trailingOnly = TRUE)[2]
freq_series <- scan(input_file,sep=',')

pdf(histo_)
hist(freq_series)
