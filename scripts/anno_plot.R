# This script takes a tables of GO terms enriched in a host group with several additional details
# and produces plots to visualize the relative enrichment of those terms
# Cyril Matthey-Doret, Laurent Casini
# 20.05.2017


in_file <- commandArgs(trailingOnly = TRUE)[1]  # Path to input file

annot_table <- read.table(in_file, header=TRUE, sep='\t')
#annot_table <- read.table('../data/annotations/HO_annot.txt',header=TRUE, sep='\t')
annot_table <- annot_table[is.finite(annot_table$oddratios),]
# Reading input file: top enriched terms in given host group

library(ggplot2)
library(gridExtra)

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
    pie(x = G$oddratios,labels = G$name,main = t_type,col=sample(seq(1,nrow(G)),size = nrow(G),replace = F) ,radius = 0.5)
  }
}


gGOpie <- function(GO_set){
  # Function is not finished. Use GOpie instead
  # This funtion takes a table of annotations as produced by
  # GO_enrich.py and produces piechart to visualize it.
  pies <- list()
  for(t_type in levels(GO_set$term_type)){
    G <- GO_set[order(GO_set$oddratios,decreasing = TRUE),]
    G <- G[G$term_type==t_type,]
    print(G)
    bp <-ggplot(data = G, aes_string(x='0',y="oddratios",fill="name")) + geom_bar(width=1,stat='identity')
    pie <- bp + coord_polar('y',start=0) + theme_minimal()
    pies[[as.character(t_type)]] <- pie
  }
  grid.arrange(grobs=pies)
}

out_file <- paste0('plots/',strsplit(basename(in_file),split = '_',fixed = TRUE)[[1]][1],'_GOpie.pdf')
pdf(out_file,width = 12,height = 18)
GOpie(annot_table)
dev.off()