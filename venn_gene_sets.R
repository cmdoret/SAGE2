setwd("~/Documents/Master/sem_2/SAGE2/src/") 
library(VennDiagram)
gene_sets <- read.csv(file = "../data/gene_sets/gene_number.csv",header=F)
B <- gene_sets$V2[gene_sets$V1=='B']+
  gene_sets$V2[gene_sets$V1=='BH']+
  gene_sets$V2[gene_sets$V1=='BO']+
  gene_sets$V2[gene_sets$V1=='BHO']

H <- gene_sets$V2[gene_sets$V1=='H']+
  gene_sets$V2[gene_sets$V1=='BH']+
  gene_sets$V2[gene_sets$V1=='HO']+
  gene_sets$V2[gene_sets$V1=='BHO']

O <- gene_sets$V2[gene_sets$V1=='O']+
  gene_sets$V2[gene_sets$V1=='BO']+
  gene_sets$V2[gene_sets$V1=='HO']+
  gene_sets$V2[gene_sets$V1=='BHO']

BH <- gene_sets$V2[gene_sets$V1=='BH']+gene_sets$V2[gene_sets$V1=='BHO']
BO <- gene_sets$V2[gene_sets$V1=='BO']+gene_sets$V2[gene_sets$V1=='BHO']
HO <- gene_sets$V2[gene_sets$V1=='HO']+gene_sets$V2[gene_sets$V1=='BHO']
BHO <- gene_sets$V2[gene_sets$V1=='BHO']

pdf("Venn_gene_sets.pdf")
grid.newpage()
draw.triple.venn(euler.d = T,scaled = T, 
                 area1 = B, 
                 area2 = H, 
                 area3 = O, 
                 n12 = BH, 
                 n23 = HO, n13 = BO, 
                 n123 = BHO, category = c("Bumblebee", "Honeybee", "Outgroup"), lty = "blank", 
                 fill = c("skyblue", "yellow", "mediumorchid"))
dev.off()
