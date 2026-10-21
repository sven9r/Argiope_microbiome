.DEFAULT_GOAL := help
.NOTPARALLEL:
.PHONY: help restore saved-figures bacteria prey history variants composition community-comparison host-phylogeny sequence-distance rarefaction figures all

RSCRIPT ?= Rscript
PYTHON ?= python3
R = $(RSCRIPT)

help:
	@echo "make restore               Restore R packages recorded in renv.lock"
	@echo "make saved-figures         Copy retained Fig1 and S2-S4 raster sources"
	@echo "make bacteria              Build bacterial state tables, S5 and S6"
	@echo "make prey                  Build prey composition tables and S9 Mantel analysis"
	@echo "make history               Build Figure3 and historical detection tables"
	@echo "make variants              Build S7 and focal 16S variant tables"
	@echo "make composition           Build Figure4; runs bacteria and prey first"
	@echo "make community-comparison  Build Figure5 with 28 displayed/25 tested profiles"
	@echo "make host-phylogeny        Draw Figure2 from the retained host-tree data"
	@echo "make sequence-distance     Build S8 from the retained sequence sources"
	@echo "make rarefaction           Build S1 using the 1500-read bacterial cutoff"
	@echo "make figures / make all    Run the complete saved-input figure workflow"
	@echo "Override runtimes with RSCRIPT=/path/to/Rscript and PYTHON=/path/to/python3."

restore:
	$(R) -e 'renv::restore(prompt = FALSE)'

# These figures have retained raster sources, not independent current renderers.
saved-figures:
	mkdir -p figures/main figures/supplement
	cp data/reference/figures/Fig01_map.png figures/main/Fig01_map.png
	cp data/reference/figures/S02_museum_success.png figures/supplement/S02_museum_success.png
	cp data/reference/figures/S03_dna_quality.png figures/supplement/S03_dna_quality.png
	cp data/reference/figures/S04_taxon_and_prey_stack.png figures/supplement/S04_taxon_and_prey_stack.png

bacteria:
	$(R) scripts/01_bacterial_profiles.R

prey:
	$(R) scripts/02_prey_profiles.R

history:
	$(R) scripts/03_historical_dusa.R

variants:
	$(R) scripts/04_symbiont_variants.R

composition: bacteria prey
	$(R) scripts/05_composition.R

community-comparison:
	$(R) scripts/06_community_comparison.R

host-phylogeny:
	$(PYTHON) scripts/07_host_phylogeny.py

sequence-distance:
	$(PYTHON) scripts/08_sequence_distance.py

rarefaction:
	$(R) scripts/09_rarefaction.R

figures: saved-figures composition history variants community-comparison host-phylogeny sequence-distance rarefaction

all: figures
