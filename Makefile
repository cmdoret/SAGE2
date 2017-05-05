include config.mk

.PHONY : all
all : Venn_gene_sets.pdf

# Generating Venn diagram
Venn_gene_sets.pdf : $(SDIR)/gene_number.csv $(SDIR)/core_number.csv scripts/venn_gene_sets.R
	Rscript scripts/venn_gene_sets.R

# Filling gene sets table
$(SDIR)/gene_number.csv : $(wildcard $(SDIR)/*_genes.txt)
	bash scripts/gene_counter.sh $@ "genes.txt"

# Filling core genes table
$(SDIR)/core_number.csv : $(wildcard $(SDIR)/*_core_set.txt)
	bash scripts/gene_counter.sh $@ "core_set.txt"

# Splitting orthologs into groups
$(SDIR)/ : $(GFAM) scripts/extract_orthol.py
	mkdir -p data/gene_sets
	python scripts/extract_orthol.py
	# Splitting orthologs into core sets of each group
	for f in B BH BHO BO H HO O;do python scripts/core_set.py $$f;done;
