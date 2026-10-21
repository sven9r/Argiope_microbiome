# Data

`Raw_data_COI_16S_2.xlsm` is the curated analysis workbook. Each marker has one row per biological specimen: 112 specimens, comprising 82 museum specimens and 30 contemporary *Argiope bruennichi*. The museum material spans 34 *Argiope* species. A1 and A2 are separate spiders.

| Worksheet | Contents | Join key |
| --- | --- | --- |
| `16S_raw_data` | Specimen metadata and 410 bacterial OTU count columns | `ID` |
| `COI_raw_data` | Specimen metadata, prey curation categories and 368 prey OTU count columns | `ID` |
| `16S_taxa` | Bacterial taxonomy annotations | `OTU` |
| `COI_taxa` | Prey taxonomy annotations | `OTU` |

Counts are nonnegative integer sequence reads. The OTU columns sum to `BacteriaReads` or `PreyReads` in every curated specimen row. OTU numbers are marker-specific: bacterial OTU2 and prey OTU2 are different features. The workbook contains 194,831 curated bacterial reads and 404,313 curated prey reads. Analyses read stored worksheet values and do not run macros.

`ID` is the biological specimen key; `MuseumID` is not necessarily unique. `Taxon` retains the source spelling of the host name. In the curated workbook, `Status` is Museum or Contemporary. In the original metadata workbook, `Age` holds that collection stratum and `Status` instead records maturity. `Year` records collection year, with unavailable values stored as `NA` or empty cells. Recorded years range from 1871 to 2021; five specimens have no year.

`TotalOTU` is the number of curated OTUs with at least 20 reads for bacteria or at least 100 reads for prey. It is not a count of OTUs after rarefaction. Original `DominantSymbiont` labels are preserved, but contemporary tests use the explicit read-composition rule in the scripts.

## Original marker tables and sequences

`raw/Tim16S.tab.txt` and `raw/TimCOI.tab.txt` are tab-separated feature-by-library count tables, with 458 bacterial and 1,122 COI features, respectively. Each has 121 library columns: 112 specimen libraries, seven c-suffix libraries and the two summaries labelled Ex and PCR. These library columns must not be counted as 121 biological specimens. Leading zeroes in raw TH labels are normalized when joining to curated IDs; the c-suffix columns remain distinct libraries.

The corresponding `*.OTu.fas` files preserve original OTU sequences. The `*.out_taxids.txt` files preserve original taxonomy exports. `raw/Argiope_Metadaten.xlsx` retains the original metadata. The raw and curated feature sets differ: a feature absent from the curated workbook is not necessarily absent from the original reads.

## Figure inputs

`reference/arg_sample_map.csv` preserves the original Arg1–Arg30 display identities. Use its original `ID` to join specimens. The display order does not change when a profile fails a read-depth cutoff.

`reference/figures/` contains the retained map and supplementary images whose original plotting workflows are not reconstructed here. `make map` runs `scripts/01_map.R` to copy the saved Figure 1 image; `make saved-figures` copies S2–S4. The original map-drawing script is unavailable. Supplementary Figure 1 has a separate current rarefaction calculation; see the root README and Makefile for the exact target.

The compact saved species summaries, pruned host tree and bacterial alignment remain at their analysis paths under `output/`. They are inputs to Figure 2 and Supplementary Figure 8. Their sources are described in [THIRD_PARTY_SOURCES.md](../THIRD_PARTY_SOURCES.md).

## Reading the proportions

Figure 2 averages per-specimen bacterial read fractions across eligible specimens; each specimen has equal weight. Detection fractions count specimens passing an operational screen. Contemporary bacterial and prey composition panels use each specimen's complete marker total. The focal-variant panel instead divides OTU2, OTU3 and OTU145 by their combined read total. These denominators answer different questions and should not be interchanged.

Raw sequencing files and verified study sequence-archive accessions are not included. This repository begins with curated count tables and retained sequence/alignment inputs. The detailed Dryad README describes control records, processing evidence, field definitions and remaining primary-record gaps.
