source("https://bioconductor.org/biocLite.R")
biocLite()

library("topGO")
library("org.Ce.eg.db")

geneList <- (rep(0, times=length(row.names(processed.table))))
names(geneList) <- row.names(processed.table)
geneList[row.names(results.de)] <- 1
geneList = as.factor(geneList)

BPdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, mapping="org.Ce.eg.db", ID="Ensembl")
MFdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.org, mapping="org.Ce.eg.db", ID="Ensembl")

resultBP <- runTest(BPdata, algorithm = "elim", statistic = "fisher")
resultMF <- runTest(MFdata, algorithm = "elim", statistic = "fisher")

GenTable(BPdata, classic = resultBP, topNodes = 100)
GenTable(MFdata, classic = resultMF, topNodes = 100)

write.table(row.names(results.de), file="foreground.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(row.names(processed.table), file="background.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

#BP stands for biological process; MF stands for molecular function; CC stands for Cellular component.
#NOTES: it is possible to change algorithm = "classic" to algorithm = "elim" in order to remove redundancy in the gene ontology graph
#NOTES: Changing topNodes to a larger values (or to length(score(resultBP)) or MF) allows to see more results. 
