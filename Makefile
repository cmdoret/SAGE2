include config.mk

.PHONY : all
all : Venn_gene_sets.pdf

# Generating Venn diagram
Venn_gene_sets.pdf : $(SDIR)/gene_number.csv venn_gene_sets.R
	Rscript venn_gene_sets.R
#
$(SDIR)/gene_number.csv : $(SDIR)/*genes.txt
	bash gene_counter.sh

$(SDIR)/*genes.txt : $(GFAM) extract_orthol.py
	mkdir -p data/gene_sets
	python extract_orthol.py

$(SDIR)/*core_set.txt : $(DIR)/*genes.txt
	for f in "B BH BHO BO";do python core_set.py;done
