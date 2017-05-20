# This script takes an ortholog table and a tab separated list of gene annotations. 
# It computes the most frequent annotation for each family and outputs a table where
# each line is a gene family with its corresponding annotation.

# Laurent Casini, Cyril Matthey-Doret
# 09.05.2017

filepath <- commandArgs(trailingOnly = TRUE)[1]  # Path to input file given as argument
################
# Loading data #
################

annotable <- read.table("./data/all_strains_annot.tsv",header=T)  #Loading annotation table
con  <- file(filepath, open = "r")  # Opening connection to input file


##############
# Processing #
##############

out_table <- data.frame(annotation=NA)  # initiating annotation table
while (length(oneLine <- readLines(con, n = 1,encoding = 'UTF-8', warn = FALSE)) > 0) {
  # Reading input file line by line
  gene_fam <- unlist(strsplit(oneLine, "\t"))  # Splitting gene family (line) into strain|gene doublets
  genfamfunc <- c()  # (Re)initializing list of functions found in the gene family
  for(d in gene_fam){
    # Iterating over strain|gene doublets
    gene <- strsplit(d, split = "|",fixed = T)[[1]][2]
    # Excluding strain name
    gene <- paste('fig',gene,sep='|')
    # Appending 'fig|' at the begining of the gene name to match annotation table format
    tmp_genfunc <- as.character(annotable$function.[annotable$feature_id==gene])
    # extracting function of corresponding gene in annotation table
    genfamfunc[tmp_genfunc] <- ifelse(test = is.numeric(genfamfunc[tmp_genfunc]),
                                      yes = genfamfunc[tmp_genfunc]+1,
                                      no = 1)
    # Incrementing number of occurence of current function in gene family by 1 (setting to 1 if never seen)
  }
  tmp_gene <- ifelse(length(genfamfunc) > 0, names(genfamfunc[which(genfamfunc==max(genfamfunc))]), NA)
  # If there was no function, returning NA, otherwise returning most frequent function
  out_table <- rbind(out_table,tmp_gene)  # Adding most frequent function in gene family to output table
}
close(con)  # Closing file connection

###########
# Writing #
###########

out_name <- strsplit(basename(filepath),".",fixed=T)[[1]][1]  # Extracting input file group (B, BH, BHO, H, HO, O)
out_table <- out_table[!is.na(out_table$annotation),]  # Removing first row containing NAs (used to initiate dataframe)
write.table(file = paste0('./data/annotations/',out_name,'_annotations.txt'),out_table,col.names = F,row.names = F,quote = F)
print(paste0("Wrote annotation file for ",out_name,". Contained ",length(out_table),
             " annotated gene families. Found ",length(unique(out_table)),  " unique functions."))
# Writing functions to file
