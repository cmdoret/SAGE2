include config.mk

.PHONY : all
all : Venn_gene_sets.pdf $(ANNOT)

# Generating Venn diagram
Venn_gene_sets.pdf : $(SDIR)/gene_number.csv $(SDIR)/core_number.csv scripts/venn_gene_sets.R
	Rscript scripts/venn_gene_sets.R

# Filling gene sets table
$(SDIR)/gene_number.csv : $(SDIR)/B_genes.txt
	bash scripts/gene_counter.sh $@ "genes.txt"

# Filling core genes table
$(SDIR)/core_number.csv : $(SDIR)/B_core_set.txt
	bash scripts/gene_counter.sh $@ "core_set.txt"

# Splitting orthologs into groups
$(SDIR)/B_genes.txt : $(GFAM) scripts/extract_orthol.py
	mkdir -p data/gene_sets
	python scripts/extract_orthol.py

# Extracting core sets (present in every strain of their group)
$(SDIR)/B_core_set.txt : $(GFAM) scripts/core_set.py
	mkdir -p data/gene_sets
	for f in B BH BHO BO H HO O;do python scripts/core_set.py $$f;done;

# Extract annotations (Long runtime!)
$(ANNOT) : scripts/ortholog_annotator.R $(SDIR)/B_genes.txt $(SDIR)/B_core_set.txt
	mkdir -p $(ANNOT)
	for f in $(SDIR)/*.txt; do Rscript $< $$f;done
