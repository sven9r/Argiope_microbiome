.DEFAULT_GOAL := help
.NOTPARALLEL:
.PHONY: help restore figure2 figure3 figure4 figure5 main-figures supplementary-figures bacteria dispersion prey history variants composition community-comparison host-phylogeny sequence-distance rarefaction processing-success dna-quality species-composition figures all

RSCRIPT ?= Rscript
PYTHON ?= python3
R = $(RSCRIPT)

help:
	@echo "Figure 1 was created in ArcGIS: figures/main/Fig01_map.png"
	@echo "make figure2                Host phylogeny and bacterial profiles"
	@echo "make figure3                Historical DUSA-like detections"
	@echo "make figure4                Bacterial and prey composition"
	@echo "make figure5                Community and prey-order comparisons"
	@echo "make main-figures           Generate Figures 2–5"
	@echo "make rarefaction            Supplementary Figure S1"
	@echo "make processing-success     Supplementary Figure S2"
	@echo "make dna-quality            Supplementary Figure S3"
	@echo "make species-composition    Supplementary Figure S4"
	@echo "make bacteria               Supplementary Figure S5"
	@echo "make dispersion             Supplementary Figure S6"
	@echo "make variants               Supplementary Figure S7"
	@echo "make sequence-distance      Supplementary Figure S8"
	@echo "make prey                   Supplementary Figure S9"
	@echo "make supplementary-figures  Generate S1–S9"
	@echo "make all                    Generate all R/Python figures"
	@echo "make restore                Restore R packages from renv.lock"

restore:
	$(R) -e 'renv::restore(prompt = FALSE)'

host-phylogeny:
	$(PYTHON) scripts/02_host_phylogeny.py

history:
	$(R) scripts/03_historical_dusa.R

composition: bacteria prey
	$(R) scripts/04_bacteria_prey_composition.R

community-comparison:
	$(R) scripts/05_community_comparison.R

figure2: host-phylogeny
figure3: history
figure4: composition
figure5: community-comparison

main-figures: figure2 figure3 figure4 figure5

rarefaction:
	$(R) scripts/S01_rarefaction.R

processing-success:
	$(R) scripts/S02_processing_success.R

dna-quality:
	$(R) scripts/S03_dna_quality.R

species-composition:
	$(R) scripts/S04_species_composition.R

bacteria:
	$(R) scripts/S05_bacterial_profiles.R

dispersion: community-comparison

variants:
	$(R) scripts/S07_symbiont_variants.R

sequence-distance:
	$(PYTHON) scripts/S08_sequence_distance.py

prey:
	$(R) scripts/S09_prey_profiles.R

supplementary-figures: rarefaction processing-success dna-quality species-composition bacteria dispersion variants sequence-distance prey

figures: main-figures supplementary-figures

all: figures
