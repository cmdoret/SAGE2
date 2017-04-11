include config.mk
 
.PHONY : all
all : Venn_gene_sets.pdf

Venn_gene_sets.pdf : $(SDIR)/gene_number.csv venn_gene_sets.R
	Rscript venn_gene_sets.R

$(SDIR)/gene_number.csv : $(SDIR)/*genes.txt
	bash gene_counter.sh

$(SDIR)/*genes.txt : $(GFAM) extract_orthol.py
	python extract_orthol.py