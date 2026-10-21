.DEFAULT_GOAL := help
.NOTPARALLEL:
.PHONY: help restore map figure1 figure2 figure3 figure4 figure5 main-figures supplementary-figures saved-figures bacteria prey history variants composition community-comparison host-phylogeny sequence-distance rarefaction figures all

RSCRIPT ?= Rscript
PYTHON ?= python3
R = $(RSCRIPT)

help:
	@echo "make restore               Restore R packages recorded in renv.lock"
	@echo "make figure1 / make map    Export the retained Figure 1 map (01_map.R)"
	@echo "make figure2               Draw the host phylogeny (02_host_phylogeny.py)"
	@echo "make figure3               Build the historical timeline (03_historical_dusa.R)"
	@echo "make figure4               Build bacterial/prey composition (04_bacteria_prey_composition.R)"
	@echo "make figure5               Build community comparisons (05_community_comparison.R)"
	@echo "make main-figures          Build Figures 1–5 in order, including their prerequisites"
	@echo "make supplementary-figures Build S1–S9; S2–S4 use retained raster sources"
	@echo "make saved-figures         Copy retained S2–S4 raster sources"
	@echo "make bacteria              Build bacterial state tables, S5 and S6"
	@echo "make prey                  Build prey composition tables and S9 Mantel analysis"
	@echo "make variants              Build S7 and focal 16S variant tables"
	@echo "make sequence-distance     Build S8 from the retained sequence sources"
	@echo "make rarefaction           Build S1 using the 1500-read bacterial cutoff"
	@echo "make figures / make all    Run main figures, then supplementary figures"
	@echo "Override runtimes with RSCRIPT=/path/to/Rscript and PYTHON=/path/to/python3."

restore:
	$(R) -e 'renv::restore(prompt = FALSE)'

map:
	$(R) scripts/01_map.R

host-phylogeny:
	$(PYTHON) scripts/02_host_phylogeny.py

history:
	$(R) scripts/03_historical_dusa.R

# Figure 4 reads the bacterial/prey exports produced by the S5/S6 and S9 scripts.
composition: bacteria prey
	$(R) scripts/04_bacteria_prey_composition.R

community-comparison:
	$(R) scripts/05_community_comparison.R

figure1: map
figure2: host-phylogeny
figure3: history
figure4: composition
figure5: community-comparison

main-figures: figure1 figure2 figure3 figure4 figure5

# S2–S4 have retained raster sources rather than current independent renderers.
saved-figures:
	mkdir -p figures/supplement
	cp data/reference/figures/S02_museum_success.png figures/supplement/S02_museum_success.png
	cp data/reference/figures/S03_dna_quality.png figures/supplement/S03_dna_quality.png
	cp data/reference/figures/S04_taxon_and_prey_stack.png figures/supplement/S04_taxon_and_prey_stack.png

rarefaction:
	$(R) scripts/S01_rarefaction.R

bacteria:
	$(R) scripts/S05_S06_bacterial_profiles.R

variants:
	$(R) scripts/S07_symbiont_variants.R

sequence-distance:
	$(PYTHON) scripts/S08_sequence_distance.py

prey:
	$(R) scripts/S09_prey_profiles.R

supplementary-figures: rarefaction saved-figures bacteria variants sequence-distance prey

figures: main-figures supplementary-figures

all: figures
