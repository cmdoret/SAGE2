# This script produces a histogram given a series of numbers. In this context, 
# each number in the series correspond to a gene family and its value represent 
# the number of strains in which it is present.
# Cyril Matthey-Doret, Laurent Casini
# 16.05.2017

input_file <- commandArgs(trailingOnly = TRUE)[1]  # Input file path as argument
out_name <- strsplit(basename(input_file),split = '_')[[1]][1]  # Extracting group
freq_series <- scan(input_file,sep=',')  # Reading input file
col_dict <- c(B='#87cdebff',
              BH='#d7ee4fff',
              BHO='#c6969aff',
              BO='#a87cdbff',
              H='#ffff00ff',
              HO='#d08e8cff',
              O='#b954d3ff')
# Colors matching each group

pdf(paste0('plots/',out_name,"_histo.pdf"))
hist(freq_series,col=col_dict[out_name],
     breaks=length(unique(freq_series)),xlab = 'Number of strains')
# Write histogram as pdf file
dev.off()
print(paste0("Generated histogram for group ", out_name))
