include config.mk

freqs=$(patsubst $(SDIR)%_genes.txt, $(FREQ)/%_gfreq.txt, $(GROUPS))

.PHONY : all
all : plots/Venn_gene_sets.pdf plots/B_GOpie.pdf plots/B_histo.pdf

# Generating Venn diagram
plots/Venn_gene_sets.pdf : scripts/venn_genfam.R $(SDIR)/gene_number.csv $(SDIR)/core_number.csv 
	mkdir -p plots
	Rscript $<

# Filling gene sets table
$(SDIR)/gene_number.csv : $(SDIR)/B_genes.txt
	bash scripts/genfam_count.sh $@ "genes.txt"

# Filling core genes table
$(SDIR)/core_number.csv : $(SDIR)/B_core_set.txt
	bash scripts/genfam_count.sh $@ "core_set.txt"

# Splitting orthologs into groups
$(SDIR)/B_genes.txt : $(GFAM) scripts/extract_orthol.py
	mkdir -p data/gene_sets
	python2 scripts/extract_orthol.py

# Extracting core sets (present in every strain of their group)
$(SDIR)/B_core_set.txt : $(GFAM) scripts/core_set.py
	mkdir -p data/gene_sets
	for f in $(GROUPS);do python2 scripts/core_set.py $$f;done;

# Calculating gene families frequencies
$(FREQ)/B_gfreq.txt : scripts/genfam_freq.py $(SDIR)/B_genes.txt
	mkdir -p $(FREQ)
	for f in $(GROUPS);do python2 $< $(SDIR)/"$$f"_genes.txt;done;

# Visualizing frequencies using histograms
plots/B_histo.pdf : scripts/histo_freq.R $(FREQ)/B_gfreq.txt
	mkdir -p plots
	for f in $(GROUPS);do Rscript $< $(FREQ)/"$$f"_gfreq.txt;done;

# Processing eggNOG annotations
$(ANNOT)/B_annot.txt : scripts/GO_enrich.py $(SDIR)/B_genes.txt
	mkdir -p $(ANNOT)
	for f in $(GROUPS);do python2 $< $(SDIR)/"$$f"_genes.txt;done;

# Generating piechart from oddratios of significant GO terms
plots/B_GOpie.pdf : scripts/anno_plot.R $(ANNOT)/B_annot.txt
	mkdir -p plots
	for f in $(GROUPS);do Rscript $< $(ANNOT)/"$$f"_annot.txt;done;
