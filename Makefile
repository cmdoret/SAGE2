include config.mk

.PHONY : all
all : Venn_gene_sets.pdf

# Generating Venn diagram
Venn_gene_sets.pdf : $(SDIR)/gene_number.csv $(SDIR)/core_number.csv venn_gene_sets.R
	Rscript venn_gene_sets.R

# Filling gene sets table
$(SDIR)/gene_number.csv : $(SDIR)/*genes.txt
	bash gene_counter.sh $* $@

# Filling core genes table
$(SDIR)/core_number.csv : $(SDIR)/*core_set.txt
	bash gene_counter.csv : $* $@

# Splitting orthologs into groups
$(SDIR)/*genes.txt : $(GFAM) extract_orthol.py
	mkdir -p data/gene_sets
	python extract_orthol.py

# Splitting orthologs into core sets of each group
$(SDIR)/*core_set.txt : $(DIR)/*genes.txt
	for f in "B BH BHO BO H HO";do python core_set.py $f;done;
