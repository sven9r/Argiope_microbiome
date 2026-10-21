# Argiope bacterial profiles and prey DNA

Bacterial 16S and prey COI analyses of 112 *Argiope* spiders: 82 museum specimens across 34 species and 30 contemporary *A. bruennichi*.

## Figures and scripts

| Figure | Content | Script or software |
| --- | --- | --- |
| **1** | Sampling map | **ArcGIS**; [map](figures/main/Fig01_map.png) |
| **2** | Host phylogeny, sampled origins, bacterial composition and detections | [02_host_phylogeny.py](scripts/02_host_phylogeny.py) |
| **3** | Historical DUSA-like read detections | [03_historical_dusa.R](scripts/03_historical_dusa.R) |
| **4** | Bacterial and prey read composition by specimen | [04_bacteria_prey_composition.R](scripts/04_bacteria_prey_composition.R) |
| **5** | Bacterial/prey NMDS, PERMANOVA and prey-order comparisons | [05_community_comparison.R](scripts/05_community_comparison.R) |
| **S1** | Rarefaction curves | [S01_rarefaction.R](scripts/S01_rarefaction.R) |
| **S2** | Marker processing outcomes | [S02_processing_success.R](scripts/S02_processing_success.R) |
| **S3** | Read recovery and OTU counts by collection period | [S03_dna_quality.R](scripts/S03_dna_quality.R) |
| **S4** | Bacterial and prey composition by host species | [S04_species_composition.R](scripts/S04_species_composition.R) |
| **S5** | DUSA-like / Mesoplasma read-count log-ratio | [S05_bacterial_profiles.R](scripts/S05_bacterial_profiles.R) |
| **S6** | Bacterial community dispersion | [05_community_comparison.R](scripts/05_community_comparison.R) |
| **S7** | Focal bacterial sequence-variant proportions | [S07_symbiont_variants.R](scripts/S07_symbiont_variants.R) |
| **S8** | Reference 16S alignment and sequence distances | [S08_sequence_distance.py](scripts/S08_sequence_distance.py) |
| **S9** | Bacterial–prey Mantel comparison | [S09_prey_profiles.R](scripts/S09_prey_profiles.R) |

[00_setup.R](scripts/00_setup.R) defines the shared data import, bacterial groups, specimen identities and plotting functions.

## Data

- [Raw_data_COI_16S_2.xlsm](data/Raw_data_COI_16S_2.xlsm): specimen metadata, curated read counts and taxonomy in `16S_raw_data`, `COI_raw_data`, `16S_taxa` and `COI_taxa`.
- [data/raw/](data/raw/): original count tables, OTU sequences, taxonomy exports and collection metadata.
- [arg_sample_map.csv](data/reference/arg_sample_map.csv): Arg1–Arg30 labels linked to specimen `ID`.
- `output/species_symbionts/` and `output/phylogeny/`: Figure 2 species summaries and published host tree; the bacterial alignment and deposited reference used for S8.

See the [data guide](data/README.md), [reference sources](THIRD_PARTY_SOURCES.md) and detailed [Dryad README](dryad_README.md) for variables, methods and provenance.

## Run the analyses

Use **R 4.5.1**, **Python 3.13** and GNU Make. Package versions are recorded in `renv.lock` and `requirements.txt`.

From the repository root:

```sh
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
make restore
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
make all PYTHON=python
```

## Makefile commands

| Command | Figures generated |
| --- | --- |
| `make figure2` | Figure 2 |
| `make figure3` | Figure 3 |
| `make figure4` | Figure 4, S5 and S9 |
| `make figure5` | Figure 5 and S6 |
| `make main-figures` | Figures 2–5 and their supporting analyses |
| `make rarefaction` | S1 |
| `make processing-success` | S2 |
| `make dna-quality` | S3 |
| `make species-composition` | S4 |
| `make bacteria` | S5 |
| `make dispersion` | S6 and Figure 5 |
| `make variants` | S7 |
| `make sequence-distance` | S8 |
| `make prey` | S9 |
| `make supplementary-figures` | S1–S9 |
| `make all` | All R/Python figures |

`make` shows the commands without running an analysis. Targets run sequentially; shared analyses run once within a single Make invocation. Figure 4's target calculates its bacterial and prey tables before drawing the panels. Outputs go to `figures/main/`, `figures/supplement/` and `output/`.

To choose a different interpreter:

```sh
make figure4 RSCRIPT=/path/to/Rscript
make sequence-distance PYTHON=/path/to/python3
```

## Analysis scope

The contemporary cutoffs are 1,500 bacterial reads and 1,693 prey reads. Figure 4 shows 29 eligible profiles per marker in 30 specimen positions; 28 specimens pass both cutoffs. Figure 5 displays those 28 and tests the 26 clear focal profiles: 16 DUSA-like and 10 Mesoplasma. TH81 and TH110 are mixed. S6 uses the same 26 paired profiles.

The profile groups compare OTU2 reads with OTU3+OTU145 reads within their combined focal total. A specimen is mixed when both focal groups contribute at least 20% or there are no focal reads. Other bacterial reads do not determine this classification. Arg1–Arg30 follow descending DUSA-like focal fraction, then descending DUSA-like reads and specimen ID.

S5 and S7 use the 28 paired specimens. S7 divides OTU2, OTU3 and OTU145 reads by their combined focal total; Figure 4 uses the complete marker total. Figure 2 applies the descriptive screen of ≥100 total bacterial reads, ≥10 focal reads and ≥1% focal reads. Figure 3 shows historical any-read detections.

The analyses start from count tables, taxonomy annotations and sequence/alignment inputs. Raw-read processing and host-tree inference are outside this code. Short-marker detections do not establish infection, strain identity or transmission. Fresh package installation has not been tested.
