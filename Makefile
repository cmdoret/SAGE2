include config.mk

.PHONY : all
all : Venn_gene_sets.pdf

# Generating Venn diagram
Venn_gene_sets.pdf : $(SDIR)/gene_number.csv $(SDIR)/core_number.csv venn_gene_sets.R
	Rscript venn_gene_sets.R

# Filling gene sets table
$(SDIR)/gene_number.csv : $(SDIR)/
	bash gene_counter.sh $@ "genes.txt"

# Filling core genes table
$(SDIR)/core_number.csv : $(SDIR)/
	bash gene_counter.sh $@ "core_set.txt"

# Splitting orthologs into groups
$(SDIR)/ : $(GFAM) extract_orthol.py
	mkdir -p data/gene_sets
	python extract_orthol.py
	# Splitting orthologs into core sets of each group
	for f in B BH BHO BO H HO O;do python core_set.py $$f;done;
