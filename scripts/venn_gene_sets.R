# This script generates a Venn diagram from a table containing the number of gene families per host group. A host group is a 
# combination of 1 to 3 of the "hosts" considered here: Bumblebee (B), Honeybee (H) and Outgroup (O). Example of host group: 
# BH for Bumblebee + Honeybee. It also takes a second table containing "core_sets". These are sets of gene families that are
# present in all strain of their group. For example, the core set BH represent gene families present in every single strain
# only found in bumblebee and honeybee, but not in outgroup.
#===========================================================
# Cyril Matthey-Doret, Laurent Casini
# 29.04.2017

###################
# Data processing #
###################

library(VennDiagram)
gene_sets <- read.csv(file = "./data/gene_sets/gene_number.csv",header=F)  # Loading list of gene sets
core_sets <- read.csv(file = "data/gene_sets/core_number.csv",header=F)  # Loading core gene sets
colnames(gene_sets) = colnames(core_sets)<- c("host_group","Ngenes")  # Naming columns for convenience
Venn_area <- c("B"=0,"H"=0,"O"=0,"BH"=0,"BO"=0,"BHO"=0,"HO"=0)  # Initializing values for the 3 Venn diagram circles
orderVenn <- c("B","BH","H","BO","BHO","HO","O")   # Categories as ordered by default in VennDiagram, will break if diagram categories change
o <- order(orderVenn)  # Vector of indexes according to orderVenn
core_sets <- core_sets[o,]  # Reordering rows of core_sets dataframe according to orderVenn
rownames(core_sets)<- NULL  # Reinitializing row names

# For each circle, add all categories containing their respective host
for(g in names(Venn_area)){
  pat=paste(rep(paste0("(.*[",g,"]{1}.*)"),times=nchar(g)),collapse = "")
  # This regex pattern will match all strings containing exactly once every character in the query,
  # no matter the order in which they appear or additional characters
  Venn_area[g] = Venn_area[g] + sum(gene_sets[grep(pattern = pat,gene_sets$host_group),"Ngenes"])
}

############
# Plotting #
############

addlab <- function(lab, x, y, offset = 0) {
  # Function for adding labels to an existing venn diagram  
  # lab: text to add on label
  # x,y: coordinates of label on plot
  # offset: vertical offset from coordinates
  # Note: a positive offset displace text downwards
  grid.text(lab, unit(as.numeric(x), "npc"), 
            unit(as.numeric(y) - offset, "npc"), 
            draw = FALSE,gp = gpar(cex=0.8,col="darkblue"))
}

pdf("Venn_gene_sets.pdf")
grid.newpage()
Venn <- draw.triple.venn(euler.d = T,scaled = T, 
                 area1 = Venn_area["B"], 
                 area2 = Venn_area["H"], 
                 area3 = Venn_area["O"], 
                 n12 = Venn_area["BH"], 
                 n23 = Venn_area["HO"], n13 = Venn_area["BO"], 
                 n123 = Venn_area["BHO"], category = c("Bumblebee", "Honeybee", "Outgroup"), lty = "blank", 
                 fill = c("skyblue", "yellow", "mediumorchid"),ind = F)

## Adding a number under each label
lbls <- gList()
o <- 1 ## counter
for(i in seq(along.with=Venn)) {
  ## Check if it is a grid.text object
  if(regexpr("text", Venn[[i]]$name) > 0 && o<=7) {
    ## Write counter value under the existing label
    lbls <- gList(lbls, addlab(core_sets[o,"Ngenes"], Venn[[i]]$x, Venn[[i]]$y, 0.03))
    ## Increase the counter
    o <- o + 1
  }
}
grid.draw(Venn)
grid.draw(lbls)
dev.off()
