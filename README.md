# Argiope bacterial profiles and prey DNA

Analysis code for bacterial 16S and prey COI data from museum and contemporary *Argiope* spiders. The curated workbook contains 112 biological specimens: 82 museum specimens across 34 species and 30 contemporary *A. bruennichi*.

## Analysis files

| File in `scripts/` | Analysis and outputs |
| --- | --- |
| [00_setup.R](scripts/00_setup.R) | Shared workbook import, bacterial grouping, specimen identities and plotting functions. |
| [01_bacterial_profiles.R](scripts/01_bacterial_profiles.R) | Contemporary bacterial read states, log-ratios (S5) and dispersion (S6). |
| [02_prey_profiles.R](scripts/02_prey_profiles.R) | Prey-order composition and bacterial/prey Mantel comparison (S9). |
| [03_historical_dusa.R](scripts/03_historical_dusa.R) | Historical DUSA-like read summaries and timeline (Figure 3). |
| [04_symbiont_variants.R](scripts/04_symbiont_variants.R) | Focal OTU sequences and within-focal read proportions (S7). |
| [05_composition.R](scripts/05_composition.R) | Aligned bacterial/prey composition panels (Figure 4). |
| [06_community_comparison.R](scripts/06_community_comparison.R) | NMDS, PERMANOVA, dispersion and prey-order comparisons (Figure 5). |
| [07_host_phylogeny.py](scripts/07_host_phylogeny.py) | Host-tree display, sampled origins, bacterial composition and detection summaries (Figure 2). |
| [08_sequence_distance.py](scripts/08_sequence_distance.py) | Reference extraction and 101-base sequence comparisons (S8). |
| [09_rarefaction.R](scripts/09_rarefaction.R) | Marker rarefaction curves and panel sample counts (S1). |


## Inputs

- [Curated workbook](data/Raw_data_COI_16S_2.xlsm): specimen metadata, read counts and taxonomy in `16S_raw_data`, `COI_raw_data`, `16S_taxa` and `COI_taxa`.
- [Original marker files](data/raw/): original count tables, OTU sequences, taxonomy exports and metadata.
- [Specimen display map](data/reference/arg_sample_map.csv): fixed Arg1–Arg30 identities, joined by biological specimen `ID`.
- `data/reference/figures/`: retained Figure 1 and S2–S4 image sources.
- Selected files under `output/species_symbionts/` and `output/phylogeny/`: saved Figure 2 summaries, published host topology, bacterial alignment and deposited reference record required by the corresponding renderers.

The [data guide](data/README.md) explains keys and count tables. [Reference sources](THIRD_PARTY_SOURCES.md) identifies the published tree and GenBank record. The detailed [Dryad README](dryad_README.md) describes the separate data package.

## Software

Use **R 4.5.1**, **Python 3.13** and GNU Make. The R lockfile records the packages used by these scripts, including `ggflags` from GitHub. Python requirements record the versions used for validation.

From the repository root:

```sh
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org")'
make restore
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
make all PYTHON=python
```

R will activate the project environment through `.Rprofile`. Restoring packages requires internet access; the analyses use the retained local inputs. Individual scripts read paths relative to the repository root.

## Makefile

`make` prints the available targets. The default does not start an analysis. `make all` and `make figures` run the same workflow; generated figures go to `figures/main/` and `figures/supplement/`, and analytical tables go to `output/`. Generated results are excluded from Git.

| Command | What it runs or writes |
| --- | --- |
| `make restore` | Restores R packages from `renv.lock`; it does not install Python packages. |
| `make saved-figures` | Copies retained Figure 1 and S2–S4 raster inputs to their publication paths. |
| `make bacteria` | Recalculates bacterial states and S5/S6. |
| `make prey` | Recalculates prey composition and S9. |
| `make composition` | Runs `bacteria`, then `prey`, then merges Figure 4 by specimen ID. |
| `make history` | Recalculates historical summaries and Figure 3. |
| `make variants` | Recalculates focal-variant summaries and S7. |
| `make community-comparison` | Runs Figure 5 with strict grouping, rarefied prey-order data, BH adjustment and starting seed 100. |
| `make host-phylogeny` | Redraws Figure 2 from the saved summaries and published host tree. |
| `make sequence-distance` | Checks the saved reference gene and alignment, exports pairwise distances and draws S8. |
| `make rarefaction` | Recalculates S1, with contemporary cutoffs of 1,500 bacterial and 1,693 prey reads. |
| `make all` | Runs the complete workflow above in the required order. |

Targets run sequentially because several share output tables. They are phony targets: invoking one runs it again, even if its outputs already exist. Repeated dependencies run once within a single Make invocation. No target deletes data or updates the lockfile.

Choose another interpreter when needed:

```sh
make composition RSCRIPT=/path/to/Rscript
make sequence-distance PYTHON=/path/to/python3
```

## Cohorts and interpretation

The contemporary bacterial cutoff is 1,500 reads and the prey cutoff is 1,693 reads. Figure 4 retains 30 specimen positions and applies each marker's cutoff separately: 29 bacterial profiles, 29 prey profiles and 28 paired profiles. The paired comparisons use those 28 specimens. Figure 5 displays all 28 and tests the 25 clear focal profiles; TH81, TH102 and TH110 are shown as mixed. Bacterial-only dispersion uses 26 clear bacterial profiles.

S5 and S7 use the paired cohort. S7 retains all 30 display positions, with the two exclusions hatched, and divides each eligible specimen's OTU2, OTU3 and OTU145 reads by their combined focal total. Figure 4 instead uses each specimen's complete marker total.

Figure 2 uses a descriptive screen: at least 100 bacterial reads, at least 10 focal reads and at least 1% focal reads. Its composition pies give eligible specimens equal weight. Figure 3 retains the historical any-read timeline selection; it is not the same screen-positive cohort.

The workflow starts from the curated workbook and saved alignment/reference inputs. It does not reconstruct raw sequencing processing, the original map, S2–S4, host-tree inference or all original species-summary calculations. `saved-figures` copies the four retained images explicitly. Randomized analyses use recorded seeds for future reruns; previously saved manuscript results can differ. A complete fresh package installation has not been tested.

Short-marker identity, candidate read detection and read-composition classes do not establish infection, bacterial strain identity or transmission. Raw sequencing accessions and complete upstream processing records remain to be supplied for a final data release.
