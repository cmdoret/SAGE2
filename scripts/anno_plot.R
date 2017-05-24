# This script takes a tables of GO terms enriched in a host group with several additional details
# and produces plots to visualize the relative enrichment of those terms
# Cyril Matthey-Doret, Laurent Casini
# 20.05.2017

library(ggplot2)
library(gridExtra)

in_file <- commandArgs(trailingOnly = TRUE)[1]  # Path to input file
in_folder <- dirname(in_file)
# Reading input file: top enriched terms in given host group

GOcbar <- function(GO_folder){
  # Producing a series of barplots from a folder as input. 
  # Each file in the folder (host grou) will be in a separate a barplot
  # Each Significantly enriched/depleted GO term will be a bar in the barplot
  
  files <- list.files(GO_folder)  # Getting all files in input folder
  clean_labels <- c(molecular_function='Molecular function', cellular_component='Cellular component', 
                    biological_process='Biological process')  # Used later to show nice plot labels
  
  for(f in files){  # Iterating over all files (besides the first one)
    tmp_table <- read.table(paste(GO_folder,f,sep='/'), header=TRUE, sep='\t')
    # Temporary table is read from file
    group_name <- strsplit(x = f,split = '_',fixed = T)[[1]][1] # Adding group column
    tmp_table <- tmp_table[is.finite(tmp_table$oddratios),]
    # Removing annotations with infinite odd ratios, caused by 0's in contingency table
    # (i.e. GO term only in one group)
    tmp_table <- tmp_table[order(tmp_table$oddratios),]  # Sorting by oddratio
    rownames(tmp_table) <- NULL  # Resetting row indices
    
    # Extracting most enriched gene for each term type (used to add labels on graph)
    # Initiating empty table to contain rows that need to be labelled
    label_table <- data.frame(term_type=numeric(0),name=numeric(0),
                              acc=numeric(0),oddratios=numeric(0),qval=numeric(0))
    for(L in levels(tmp_table$term_type)){  # Iterating over term types
      term_rows <- tmp_table[tmp_table$term_type==L,]  # Subsetting table
      term_rows$name <- as.character(term_rows$name)
      term_rows$name[term_rows$oddratios<max(term_rows$oddratios)] <- ''
      label_table<- rbind(label_table, term_rows[,c('term_type','name',
                                                    'acc','oddratios','qval')])
      # For each term type, appending row(s) with highest enrichment value
    }
    bp <- ggplot(data=tmp_table, aes(x=reorder(acc, oddratios),y=log10(oddratios),fill=qval)) +
      geom_bar(stat="identity")
    # Using ggplot to produce barplots facetted by group. Barplots measure odd ratios of different annotations
    bp <- bp + theme_minimal() + theme(axis.text.x = element_blank()) + 
      facet_grid(~term_type, scales='free', drop=T, labeller = as_labeller(clean_labels))
    bp <- bp + geom_hline(aes(yintercept=0)) + ggtitle(group_name) +
      geom_text(data=label_table,aes(label=name, x=reorder(acc, oddratios),
               y=1.1*log10(oddratios))) + xlab("GO terms") + ylab("Odds ratio [log10]")
    # Adding a line to mark 1 (below-> depleted, above -> enriched)
    
    # Sending plot to output file
    out_file <- paste0('plots/',strsplit(basename(f),split = '_',fixed = TRUE)[[1]][1],'_GObar.pdf')
    pdf(out_file)
    print(bp)
    dev.off()
  }
}

GOcbar(in_folder)


#===========================================
# THESE FUNCTIONS ARE NOT USED ANYMORE, USE GOcbar INSTEAD

#annot_table <- read.table(in_file, header=TRUE, sep='\t')
#annot_table <- read.table('../data/annotations/HO_annot.txt',header=TRUE, sep='\t')
#annot_table <- annot_table[is.finite(annot_table$oddratios),]
#annot_table <- annot_table[annot_table$oddratios>1,]

# out_file <- paste0('plots/',strsplit(basename(in_file),split = '_',fixed = TRUE)[[1]][1],'_GOpie.pdf')
# pdf(out_file,width = 12,height = 18)
# GOpie(annot_table)
# dev.off()

GOpie <- function(GO_set){
  # This funtion takes a table of annotations as produced by
  # GO_enrich.py and produces piecharts to visualize it.
  par(mfrow=c(3,1))
  for(t_type in levels(GO_set$term_type)){
    G <- GO_set[order(GO_set$oddratios,decreasing = TRUE),]
    G <- G[G$term_type==t_type,]
    thresh_enrich <- sort(G$oddratios,decreasing = T)[10]
    G$name <- as.character(G$name)
    G$name[G$oddratios<thresh_enrich] <- ''
    pie(x = G$oddratios,labels = G$name,main = t_type,col=
          sample(seq(1,nrow(G)),size = nrow(G),replace = F) ,radius = 0.5)
  }
}


GObar <- function(GO_folder){
  # Producing a series of barplots from a folder as input. 
  # Each file in the folder (host grou) will be in a separate a barplot
  # Each Significantly enriched/depleted GO term will be a bar in the barplot
  files <- list.files(GO_folder)  # Getting all files in input folder
  full_table <- read.table(paste(GO_folder,files[1],sep='/'), header=TRUE, sep='\t')
  # Initializing table with first file
  full_table$group <- strsplit(x = files[1],split = '_',fixed = T)[[1]][1]
  # Adding group column. info extracted from filename! :(
  for(f in files[2:length(files)]){  # Iterating over all files (besides the first one)
    tmp_table <- read.table(paste(GO_folder,f,sep='/'), header=TRUE, sep='\t')
    # Temporary table is read from file
    tmp_table$group <- strsplit(x = f,split = '_',fixed = T)[[1]][1] # Adding group column
    full_table <- rbind(full_table,tmp_table)  # Temporary table appended to full table
  }
  full_table <- full_table[is.finite(full_table$oddratios),]
  # Removing annotations with infinite odd ratios, caused by 0's in contingency table
  # (i.e. GO term only in one group)
  rownames(full_table) <- NULL  # Resetting row indices
  full_table$group <- as.factor(full_table$group)  # group column as factor for facetting in ggplot
  print(full_table)
  bp <- ggplot(data=full_table, aes(x=name,y=oddratios))+geom_bar(stat="identity")
  # Using ggplot to produce barplots facetted by group. Barplots measure odd ratios of different annotations
  bp <- bp + geom_hline(aes(yintercept=1)) + facet_grid(group~term_type)
  # Adding a line to mark 1 (below-> depleted, above -> enriched)
  bp
}
