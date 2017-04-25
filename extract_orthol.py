
# The purpose of this script is to extract a list of genes that are orthologous  between bumblebee (B)- and honeybee (H)-specific
# bacterial species. It takes an ortholog table as input and output a similar table containing only B-H orthologs.
# Laurent Casini, Cyril Matthey-Doret
# 04.04.2017

#==========================================
## Loading data

# List of genomes IDs
ID_H = ['JF72','F259','JG30','JF76','F260','F261','F262','F263',
        'JF73','JF74','JF75','WB8','WB10','L185','L186',
        'L184','L183'] # Honeybee-bacterial genomes
ID_B = ['F225','F230','F233','F234','F236','F237',
        'F228','F245','F246','F247'] # Bumblebee-bacterial genomes
ID_O = ['LA14','LA2','LDB','LGAS','LHV','LJP','WANG','JG29'] # Outgroup bacterial genomes

# Loading ortholog table and creating output files
ortho_tab = open("./data/mclOutput",'r')
B_genes = open("./data/gene_sets/B_genes.txt",'w')
H_genes = open("./data/gene_sets/H_genes.txt",'w')
HO_genes = open("./data/gene_sets/HO_genes.txt",'w')
BH_genes = open("./data/gene_sets/BH_genes.txt",'w')
BO_genes = open("./data/gene_sets/BO_genes.txt",'w')
BHO_genes = open("./data/gene_sets/BHO_genes.txt",'w')
O_genes = open("./data/gene_sets/O_genes.txt",'w')
#==========================================

for line in ortho_tab:  # Each line is a gene family
    tmp = line.split("\t")  # split line by word (tab separated)
    B = False; H = False; O = False  # Is this family present in bumble/honeybees
    for word in tmp: # For each word (a word is: "strain|gene")
        gene = word.split("|")  # sparates strain from gene in a list

        if gene[0] in ID_H: # if the strain is in honeybee
            H = True
        if gene[0] in ID_B: # if the strain is in bumblebee
            B = True
        if gene[0] in ID_O: # if the strain is in outgroup
            O = True

    if B and H and O:  # if gene family in all strains
        BHO_genes.write(line)  # writing family to all file
    elif B and H and not O:
        BH_genes.write(line)
    elif B and not H and O:
        BO_genes.write(line)
    elif B and not H and not O:
        B_genes.write(line)
    elif not B and H and O:
        HO_genes.write(line)
    elif not B and H and not O:
        H_genes.write(line)
    elif not B and not H and O:
        O_genes.write(line)



ortho_tab.close();BHO_genes.close()
B_genes.close();BH_genes.close();BO_genes.close()
H_genes.close();HO_genes.close();O_genes.close()
