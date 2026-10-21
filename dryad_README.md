# Bacterial 16S and prey COI data from museum and contemporary Argiope

Data documentation prepared **21 October 2026**. This is a local deposit draft; it has not been submitted to Dryad or assigned a dataset DOI. Dataset creators, the final publication citation and project data licence must be supplied with the deposit metadata.

## Dataset description

The dataset contains short-marker bacterial and prey sequences, specimen-level read counts and collection metadata for 112 *Argiope* spiders. Museum material provides comparative host and historical coverage; contemporary *A. bruennichi* provides paired bacterial and prey profiles from one German collection. The focal bacterial features are DUSA-like OTU2 and the source-annotated *Mesoplasma* OTU3 and OTU145.

Read counts, observed detections, descriptive screening categories and predominant bacterial profiles are separate measurements. A short sequence detection does not establish tissue infection, bacterial function, transmission or absolute bacterial abundance. Prey COI reads do not directly measure prey biomass or the bacteria carried by prey.

The two markers have separate OTU namespaces. Bacterial OTU2 and prey OTU2 are unrelated identifiers. Original OTU identifiers are retained throughout; these files do not introduce a new ASV inference.

| Dataset component | Records or features |
| --- | ---: |
| Biological specimens | 112 |
| Museum specimens | 82 |
| Contemporary specimens | 30 |
| Sampled host species | 34 |
| Museum *A. bruennichi* | 12 |
| Contemporary *A. bruennichi* | 30 |
| Specimens with numeric collection years | 107 |
| Recorded collection-year range | 1871–2021 |
| Curated bacterial OTUs | 410 |
| Curated prey COI OTUs | 368 |
| Sum of curated bacterial specimen reads | 194,831 |
| Sum of curated prey specimen reads | 404,313 |
| Original bacterial OTU sequences/count-table rows | 458 |
| Original prey COI OTU sequences/count-table rows | 1,122 |
| Original count-table library columns per marker | 121 |

The 121 library columns comprise 112 specimen-labelled columns, seven additional c-suffix columns and two pooled control columns. The latter are not additional biological specimens. A1 and A2 are two separate spiders, both *A. argentata*; their shared museum catalogue number does not identify a duplicate specimen.

## Study material and geographic scope

The working manuscript identifies the museum source as the arachnological collection of the Senckenberg Research Institute, Frankfurt am Main, Germany. Contemporary female *A. bruennichi* were collected in 2021 from a meadow near Schwanheimer Dünen, Frankfurt am Main, and preserved in ethanol. Museum specimens represent multiple species, countries and collection years.

Use each specimen's recorded `Country`, `Location` and `Year` for its source context. Country and locality text are preserved as recorded, including spelling variants and territory names. These are sampled origins, not complete host distributions. Locality precision and coordinate uncertainty are not encoded in the curated count tables. A sampling map cannot supply missing specimen coordinates.

The metadata include adult, subadult and juvenile annotations. Recorded sex is `female` or `NA`; missing sex annotations should remain missing even though the manuscript describes female sampling. Five specimens have unknown collection years: A1, A2, TH52, TH79 and TH80.

## File organization

All paths below are relative to the deposit package root. `file_manifest.csv` is the complete inventory, with byte size, SHA-256, original project path, source checksum and export operation. It also identifies the worksheet behind each worksheet export.

| File or directory | Contents and intended use |
| --- | --- |
| `curated/` | Current specimen-by-OTU matrices and marker-specific taxonomy tables; the starting point for most reuse. |
| `metadata/` | Original specimen collection metadata in CSV form. |
| `raw_counts/` | Original, unfiltered OTU-by-library matrices, library identities and observed control counts. |
| `source_workbooks/` | Unaltered source workbooks, including nonrectangular sheets and original formula/macro containers. |
| `sequences/` | Original marker sequences, retained bacterial alignment/tree and complete five-sequence S8 comparison. |
| `figure_support/current/` | Selected current-run support with an identified source, including Figure 5 and S1 sidecars. |
| `figure_support/species_symbionts/` | Existing candidate-screen, species-summary, taxonomy, host-mapping and composition support. Some composition/reconciliation files retain older analysis thresholds; see the figure map. |
| `figure_support/legacy/` | Older generic tables and saved statistical outputs. They are retained as provenance, not relabelled as the current selected run. |
| `figures/main/` | Five manuscript main-figure PNGs. |
| `figures/supplement/` | Nine supplementary-figure PNGs. |
| `third_party/` | Published host-tree source and deposited DUSA reference, with provenance. |
| `legacy_context/` | Exports from the earlier curated workbook, retained for comparison only. |
| `provenance/` | Source metadata and methods records needed to interpret the deposited data. |
| `column_dictionary.csv` | File-specific column definitions, observed storage types, missing tokens and explicitly unresolved definitions. |
| `figure_file_map.csv` | Figure filename, associated support files, original source and version limitations. |
| `package_checks.json` | Export dimensions, marker row-sum checks and package scope. |
| `METHODS_AND_GAPS.md` | Detailed methods evidence and limits of reproducibility. |
| `THIRD_PARTY_SOURCES.md` | Attribution, accessions and recorded terms for reused material. |

Analysis code is maintained separately from the data deposit. The package contains neither installed software libraries nor repository inventory, cleanup, document-synchronization or packaging helpers.

## Current curated count matrices

### Bacterial 16S

`curated/16s_specimen_counts.csv` contains 112 rows and 426 columns: 16 specimen/summary fields followed by 410 OTU count columns. It is exported from `16S_raw_data` in `source_workbooks/Raw_data_COI_16S_2.xlsm`.

Each row represents one biological specimen. Each OTU cell is a nonnegative integer count of reads assigned to that bacterial feature. `BacteriaReads` equals the sum of all 410 OTU cells for that specimen. The matrix retains low-depth and zero-read specimens; inclusion in a particular analysis must be determined separately.

`curated/16s_taxonomy.csv` contains 410 rows, one per curated bacterial OTU, with these columns:

`OTU`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`.

These are retained source annotations. The source label `Unknown symbiont` corresponds to focal DUSA-like OTU2 in the analysis; it is not a newly described taxon. Unresolved ranks remain unresolved.

### Prey COI

`curated/coi_specimen_counts.csv` contains the same 112 specimen IDs, with 390 columns: 22 specimen/curation/summary fields and 368 OTU columns. It is exported from `COI_raw_data` in the same current workbook.

`PreyReads` equals the sum of the 368 curated prey OTU cells for that specimen. The matrix includes unsuccessful, ambiguous and invalid source processing categories as well as valid profiles. A row's presence does not imply that it passed a particular prey analysis.

`curated/coi_taxonomy.csv` contains 368 rows with these columns:

`OTU`, `Match`, `Seq`, `Kingdom`, `Phylum`, `Class`, `Subclass`, `Order`, `Superfamily`, `Family`, `Genus`, `Species`.

`Match` and `Seq` retain their stored values. Their numerical scaling and exact original definitions have not been recovered. A stored `Match` value such as `100000` must not be read as a documented percentage identity. The taxonomy table does not supply a complete BLAST alignment, database date, accession trail or uncertainty model for each prey assignment.

## Specimen and curation fields

The following definitions apply to the **current curated specimen tables**. Some headers in the original metadata have different meanings, as described below.

| Field | Definition, units and interpretation |
| --- | --- |
| `Number` | Source specimen ordering number; dimensionless. Not a read count or stable cross-file join key. |
| `ID` | Unique specimen identifier: A1, A2, TH1–TH110. Primary join key between markers. |
| `MuseumID` | Museum catalogue number. Multiple specimens can share it; it is not a unique specimen key. |
| `Taxon` | Recorded spider host identification. Original spacing and abbreviations are retained, including `A.modesta`. |
| `Status` | Collection stratum: `Museum` or `Contemporary`. |
| `Maturity` | Recorded maturity: `adult`, `subadult` or `juvenil`; source spelling is retained. |
| `Sex` | Recorded sex, or `NA` when unavailable. |
| `Continent` | Recorded continent text. |
| `Country` | Recorded country/territory text. |
| `Location` | Recorded locality text. No coordinate units apply. |
| `ConditionOpisthosoma` | Recorded condition: `intact`, `damaged`, `separated`, or `NA`. |
| `Year` | Collection year, in calendar years, or source `NA`. It is not specimen age at extraction. |
| `DominantSymbiont` | Original source annotation: `Unknown symbiont`, `Mesoplasma` or `NA`. Distinct from the later explicit profile-group rule. |
| `UnknownSymbiont` | Original `dominant` / `not dominant` annotation, or `NA`. |
| `TotalOTU` in bacterial table | Number of curated bacterial OTUs with **at least 20 reads** in that unrarefied specimen profile. |
| `TotalOTU` in prey table | Number of curated prey OTUs with **at least 100 reads** in that unrarefied specimen profile. |
| `BacteriaReads` | Complete curated bacterial read total; integer reads. |
| `PreyReads` | Complete curated prey read total; integer reads. |
| `OTU1`, `OTU2`, … | Marker-specific feature read counts. The actual feature set and column order are recorded in each file and its dictionary. |
| `Category` | Original prey processing category: `valid`, `ambiguous`, `invalid` or `inadmissible`. Detailed classification criteria are incompletely documented. |
| `Green`, `Yellow`, `Orange`, `Red`, `Blue` | Stored prey curation fields. Exact scoring/category rules have not been recovered; do not assign new biological meanings to the colors. |

Both `TotalOTU` definitions agree with all 112 rows of their respective curated matrices. They are neither ordinary nonzero-feature richness nor richness recalculated after rarefaction. The prey workbook additionally retains a `COUNTIF(..., ">=100")` formula supporting that definition. The bacterial definition is established by exact correspondence of all stored values to counts at the ≥20-read threshold.

Across all 112 prey records, `Category` contains 48 `valid`, eight `ambiguous`, three `invalid` and 53 `inadmissible` values. The categories are source processing outcomes; they must not be replaced by a read-depth flag. Contemporary records comprise 29 valid and one inadmissible profile.

## Original collection metadata

`metadata/original_metadata.csv` exports the 112 rows in the `Metadata` worksheet of `source_workbooks/Argiope_Metadaten.xlsx`.

| Original header | Meaning and cross-reference |
| --- | --- |
| `Age` | `Museum` / `Contemporary` collection stratum; corresponds to curated `Status`. It is not numeric age. |
| `Sample ID` | Specimen identifier; join to curated `ID`. |
| `Museum ID` | Catalogue number; corresponds to curated `MuseumID`. |
| `Taxon` | Recorded host species. |
| `Status` | Maturity annotation; corresponds to curated `Maturity`, not curated `Status`. |
| `Sex` | Recorded sex. |
| `Continent`, `Country`, `Location` | Original recorded origin text. |
| `Year of catch` | Collection year; corresponds to curated `Year`. |

The workbook also retains `Tabelle Methodenteil`, a secondary source worksheet. Its original layout is preserved in the workbook; the rectangular metadata CSV does not replace every source worksheet.

## Original counts and library identities

`raw_counts/16s_original_counts.csv` and `raw_counts/coi_original_counts.csv` are transposed relative to the curated matrices: one row is one original OTU and each remaining column is a library. The first header is `#OTU ID`. The bacterial matrix has 458 feature rows; COI has 1,122. Both contain the same 121 library labels.

These matrices retain the larger original feature universe. Their totals are not interchangeable with `BacteriaReads` or `PreyReads`, which sum only curated features. Both original marker FASTAs are included so that feature IDs remain connected to their measured sequences.

`raw_counts/library_identity_map.csv` records marker, original library label, normalized specimen ID, library role and the scope of design verification. The additional c-suffix labels are TH10c, TH18c, TH36c, TH40c, TH50c, TH58c and TH65c. Source metadata map them to existing specimens. Their experimental level—PCR, extraction, tissue or another repeat—has not been established from laboratory records.

For all seven mapped specimens and both markers, curated retained-feature counts equal the original **base TH library**, rather than the sum of base and c-suffix libraries. Do not add c-suffix columns to the curated profiles when reproducing the present analysis.

### Control columns

The original text matrices contain `Ex` and `PCR` columns. Worksheet evidence establishes that they pool separately named source libraries: `Ex_c`, `Ex_c2`, `PCR_c`, `PCR_c2` and `PCR_c3`. `Ex` and `PCR` are pooled count columns, not two independently established physical controls.

| Marker / feature universe | Ex reads | PCR reads |
| --- | ---: | ---: |
| Bacterial, all 458 original OTUs | 5,324 | 4,181 |
| COI, all 1,122 original OTUs | 610 | 1,710 |

For every bacterial OTU, `Ex = Ex_c + Ex_c2` and `PCR = PCR_c + PCR_c2 + PCR_c3`. The five full-universe bacterial totals are 1,792; 3,532; 672; 1,481; and 2,028 reads, respectively. `raw_counts/16s_control_library_counts.csv` documents the observed separate-library counts, including focal OTUs.

Bacterial OTU2 has pooled Ex/PCR counts of 2/0. OTU3 has 2/2 and OTU145 has 0/0. These are observations, not a contamination correction. The descriptive candidate screen does not automatically subtract these counts from specimens.

The COI workbook's detailed five-library columns cover 940 of the 1,122 original OTUs. Within those 940 features, the pooled equations agree exactly and give 206 Ex and 673 PCR reads. The remaining 182 features account for 404 Ex and 1,037 PCR reads; their allocation among the five separate libraries cannot be inspected in that sheet.

The working methods report extraction and PCR negatives. Physical control identities, negative/positive roles, sample-batch coverage, additional positive controls and the original contaminant-removal algorithm still require laboratory records. The pooled columns do not certify completeness of the control design.

## Joining specimens, features and figures

1. Join the curated marker tables on `ID`, checking uniqueness on each side. They contain identical sets of 112 specimen IDs. Preserve a marker label whenever count columns are reshaped.
2. Join original metadata using `Sample ID` = curated `ID`. Rename original `Age` and `Status` before combining them with curated `Status` and `Maturity`.
3. Within each marker, join OTU counts to its own taxonomy table by `OTU`. Do not join bacterial taxonomy to prey counts with the same identifier.
4. Original raw labels such as `TH01` normalize to curated `TH1`. A terminal `c` can be removed for specimen identity lookup only. Keep the original library label and role; normalization is not permission to merge counts.
5. FASTA headers use source spelling such as `OTu2`; normalize case only when matching to `OTU2`. Retain a record of normalization. Do not truncate or renumber feature IDs.
6. Join plotted Arg labels through the specimen-ID map associated with that figure. `Arg1` is a display label, not a universal specimen ID. The Figure 4 union retains Arg1–Arg30 positions; optional matched companions can use different sequential labels.
7. For distance-based analyses, restrict both markers to the same eligible specimen IDs and put matrices in exactly the same row order before comparing them.

A1/A2 remain specimen records in these joins. `MuseumID`, row order, country, host species and the displayed Arg label do not replace `ID` as the specimen key.

## Analysis inclusion and denominators

The same curated profiles support several questions. There is no single universal read cutoff.

| Scope | Bacterial minimum | Prey minimum | Eligible specimens |
| --- | ---: | ---: | --- |
| Descriptive host/candidate screen | 100 | Not required | 100: 70 museum + 30 contemporary |
| Global bacterial recovery/composition | 507 | Not required | 56: 27 museum + 29 contemporary |
| Global prey recovery | Not required | 1,133 | 59: 30 museum + 29 contemporary |
| Global prey species stacks | Not required | 1,133 and `Category == "valid"` | 48: 19 museum + 29 contemporary |
| Contemporary bacterial profile | 1,500 | Not required | 29 of 30 |
| Contemporary prey profile | Not required | 1,693 | 29 of 30 |
| Paired contemporary comparison | 1,500 | 1,693 | 28 of 30 |
| Clear two-state paired comparison | Same paired depths and explicit profile rule | Same | 25: 15 DUSA-like + 10 *Mesoplasma* |

TH84 has 4,375 bacterial and 746 prey reads and fails the paired prey cutoff. TH99 has 182 bacterial and 3,207 prey reads and fails the paired bacterial cutoff. TH102 has 1,861 bacterial and 4,061 prey reads and is retained at 1,500 bacterial reads. Older 2,506-read analyses excluded TH102 and used a 27-specimen paired cohort; their outputs are not the current 28-specimen analysis.

Figure 4 retains 30 specimen positions and marker-specific eligibility, giving 29 eligible profiles per panel and 28 with both markers. Figure 5, S5 and the current paired S7 display use the 28-specimen intersection; S7 retains 30 fixed display positions. Failed-marker or paired-exclusion hatches are missing/excluded information, not zero proportions.

### Candidate detection screen

A profile is eligible for descriptive screening when its complete curated bacterial total is ≥100 reads. A candidate group passes when it also has **at least ten candidate reads and at least 1% of the complete curated bacterial total**. `Trace-only` means candidate reads are present in an eligible profile but one or both candidate thresholds fail. `No candidate reads` means zero recorded reads in an eligible profile. Below-eligibility profiles are separate.

Group membership is explicit in `candidate_group_definitions.csv`; DUSA-like is OTU2 and *Mesoplasma* is OTU3 + OTU145. OTU287 remains a separate close sequence relative. Named, `_hint`, `_unresolved` and unnamed-Candidatus categories preserve different evidence levels. The catalogue is an operational annotation of retained labels, not a new comprehensive database search.

Candidate screens are independent: one specimen can pass both the DUSA-like and *Mesoplasma* screens. These descriptive fractions are not validated infection prevalence or mutually exclusive dominance classes.

### Relative profile groups and log-ratio

The contemporary profile rule flags a specimen as mixed/ambiguous when any of these holds:

- Its largest individual bacterial OTU is outside OTU2, OTU3 and OTU145.
- DUSA-like and combined *Mesoplasma* each contribute at least 20% of complete bacterial reads.
- `abs(log10((OTU2 + 1) / (OTU3 + OTU145 + 1))) < 0.5`.

For unambiguous profiles, the group follows the focal group containing the largest individual OTU. The one-read pseudocount is added to both group counts before the log-ratio. The rule uses the largest **individual OTU**, rather than simply comparing the two summed focal groups. The paired cohort contains 15 DUSA-like, ten *Mesoplasma* and three mixed profiles: TH81, TH102 and TH110.

The bacterial-only set contains 26 clear profiles, including TH84. A saved dispersion calculation based on 25 clear **paired** profiles is a different cohort and must not be described as the 26-profile bacterial-only calculation.

### Percentage denominators

| Quantity | Denominator and scale |
| --- | --- |
| Per-specimen bacterial composition | All 410 curated bacterial features in that specimen; fraction 0–1 or percent 0–100 as named. |
| Per-specimen prey composition | All 368 curated prey features in that specimen; percent 0–100 in detailed composition tables. |
| Figure 2 composition | Arithmetic mean of eligible per-specimen bacterial fractions within species × collection stratum; equal specimen weight. |
| Candidate screen fraction | Candidate-group reads / complete curated bacterial total for that specimen. |
| Candidate screen-positive proportion | Number passing / number eligible within the stated species × collection stratum. |
| Focal-variant composition, S7 | OTU-specific reads / `(OTU2 + OTU3 + OTU145)` for that specimen. Other bacteria are excluded. |
| Within-group variant fraction | OTU reads / the named focal group's reads; undefined when that group has no reads. |
| Global species stacks, S4 | Counts rarefied to 507 bacterial or 1,133 prey reads, pooled within host species and normalized to 100%. |

Do not substitute pooled-read fractions for arithmetic means of specimen fractions. A focal-only S7 percentage is not the fraction of the complete bacterial community.

## Figure support and saved analyses

`figure_file_map.csv` associates each image with support and source. Historical filenames such as `S04_logratio_*` or `S07_detailed_diet_*` can refer to analyses now displayed as S5 or Figure 4B; use the map rather than inferring figure numbering from filenames.

| Figure | Data represented and interpretation |
| --- | --- |
| Figure 1 | Sample origins of museum and contemporary material. Overlapping symbols are not distinct localities or complete ranges. |
| Figure 2 | Published host topology; separate museum/contemporary *A. bruennichi* rows; sampled origins; ≥100-read bacterial composition means and candidate screens. |
| Figure 3 | Retained historical OTU2 any-read observations. The newer screen-positive timeline is a separate selection and is not silently substituted for the retained image. |
| Figure 4 | Contemporary bacterial/prey composition aligned by specimen; complete-marker denominators, 29 profiles per marker and 30 fixed positions. |
| Figure 5 | Paired ordinations and prey-order comparisons; 28 displayed profiles and 25 clear profiles for two-state tests in the 1,500-read run. |
| S1 | Rarefaction context; current 1,500/1,693 contemporary depths and global 507/1,133 depths. |
| S2 | Original processing outcomes, distinct from paired-profile eligibility. |
| S3 | Read recovery and thresholded source `TotalOTU` counts by collection period. |
| S4 | Species-pooled rarefied composition using global depths; museum/contemporary material pooled within species. |
| S5 | Focal bacterial log-ratio with a +1 read pseudocount, 28 paired profiles and explicit mixed profiles. |
| S6 | Bacterial-only dispersion display; 26 clear profiles at the current bacterial depth. Figure 5's paired dispersion sidecar contains 25 and is not its source. |
| S7 | Three original focal OTU variants divided by combined focal reads; 28 paired profiles in 30 positions. The filename's “haplotypes” does not denote host haplotypes. |
| S8 | Ten unordered sequence comparisons among five 101-base fragments. |
| S9 | Paired bacterial–prey distance comparison; specimen-label permutations, not independent distance-pair replicates. |

Current selected sidecars take precedence over older generic outputs. Older files can retain different thresholds, transformations, groups and stochastic runs. A model table belongs to a PNG only when the selected source/run supports the link. The file map records missing or superseded sidecars rather than attaching an old result to a new image.

### Collection periods and source OTU counts

S3 uses right-closed intervals: year ≤1900; 1900 < year ≤1950; 1950 < year ≤2000; 2000 < year ≤2020; year >2020. All 30 records in the last interval are contemporary specimens from 2021. No 1900 or 2020 observation is present. The retained image's first `<1900` label is narrower than the coded interval but changes no current specimen membership.

The read panels include 107 dated specimens. The bacterial-OTU panel includes 55 dated profiles with ≥507 bacterial reads; the prey-OTU panel includes 58 with ≥1,133 prey reads. Period sample counts are 10/20/29/18/30 for reads, 1/3/12/10/29 for bacterial OTU counts, and 2/2/10/15/29 for prey OTU counts.

Original Kruskal–Wallis annotations and Dunn/Bonferroni letter relations were checked against these source variables. They concern thresholded unrarefied `TotalOTU`, not recalculated rarefied richness. Historical software versions and literal compact-letter assignment order have not been fully recovered.

S4 includes 56 bacterial profiles across 19 species and 48 valid prey profiles across 12 species. Its *A. bruennichi* labels are **34 bacterial and 35 prey profiles**. An earlier local note read the bacterial label as 44; that reading was corrected. It is not evidence of a 44-profile dataset.

### Community distances and tests

Current contemporary analyses average Bray–Curtis dissimilarities across 100 subsamples to 1,500 bacterial or 1,693 prey reads. These distances are separate from unrarefied proportions. A paired Mantel run contains 28 specimens, generating 378 unordered distances. Those 378 pairs are derived comparisons from 28 specimens, not 378 biological replicates.

The 21 October 1,500-read S9 run uses a Spearman Mantel statistic and 999 specimen-label permutations for positive association. Its retained run records seeds 100, 101 and 102 for bacterial distances, prey distances and Mantel, respectively. The plotted linear fit is descriptive, not the Mantel statistic. Consult the figure map for whether complete matching current-image sidecars are included; the legacy 27-specimen table is not this run.

Figure 5 uses 25 clear paired profiles for the prey comparison, 999 PERMANOVA permutations and seven unpaired prey-order tests with Benjamini–Hochberg adjustment. Bacterial features define the profile groups; a bacterial PERMANOVA separation is not independent evidence of a group difference. Saved scores/statistics retain their selected run rather than certifying all older tables.

## Sequence files and bacterial reference comparison

| File in `sequences/` | Sequence content |
| --- | --- |
| `original_16s_otus.fasta` | All 458 original bacterial representatives. All except OTU116 are 87–111 bases; OTU116 is 211 bases. |
| `original_coi_otus.fasta` | All 1,122 original prey COI sequences, each 71 bases. |
| `retained_16s_no_OTU116.fasta` | 457 bacterial sequences after the single stated exclusion of OTU116. |
| `retained_16s_forced_alignment.fasta` | The same 457 sequences in a 156-column reference-guided alignment. Gaps are alignment characters, not measured nucleotides. |
| `retained_16s_fasttree_unrooted.nwk` | Unrooted bacterial tree from 102 alignment columns with ≥80% non-gap coverage. |
| `S08_measured_101bp_fragments.fasta` | Five complete 101-base fragments: reference, OTU2, OTU287, OTU3 and OTU145. |

OTU116 is absent from the retained alignment/tree and current curated matrix but preserved in original evidence. Its unusual length and source Streptophyta hint motivated exclusion. Of the remaining 457 sequences, 47 lack a curated count column; that absence does not mean zero original reads. Twelve shared TH2 count cells are zero in the curated workbook and nonzero in the raw table; the difference table preserves them without inventing their curation rationale.

The DUSA reference is [CP070274.1](https://www.ncbi.nlm.nih.gov/nuccore/CP070274.1). Its annotated 16S gene is extracted in gene orientation from `complement(553873..555402)` of the genome and is 1,530 bases long. OTU2 uniquely matches gene positions 1060–1160, inclusive. Genome and oriented-gene coordinates are different coordinate systems.

S8 distances are uncorrected mismatches divided by 101 homologous A/C/G/T sites. OTU2 matches 101/101 reference bases; OTU287 differs at four sites. OTU3/OTU145 differ from OTU2 at 19/15 sites and from each other at four sites. No species cutoff is imposed. Exact fragment identity does not establish identical genomes, strains or infection.

The 457-tip bacterial tree is distinct from the published host tree. Short amplicons support exploratory sequence relationships; deep relationships and species identity should not be inferred from this tree alone.

## Molecular methods and evidence limits

The working methods describe gut-tissue removal after 0.15% sodium hypochlorite treatment, an ethanol rinse and drying; Qiagen Puregene extraction with overnight Proteinase K digestion and GlycoBlue; short-marker PCRs; dual indexing in a second PCR; AMPure XP cleanup; and Illumina MiSeq 300-cycle sequencing. These are reported methods. Counts and FASTAs do not independently certify each laboratory step.

Reported primer pairs are 16S_Museum_F1 / 16S_Museum_R1 and NoSpid_edit_F1 / Laurelin_edit_R1. The latter correspond to Spidprey_F / Spidprey_R in Melcher et al. (2024), [doi:10.1111/een.13382](https://doi.org/10.1111/een.13382). COI sequences are 71 bases after primer removal. The original 16S table reports 122 bp, whereas the verified focal retained fragment is 101 bases. Original primer-design provenance and the length definition remain unresolved; neither value should be silently substituted for the other.

The reported processing lineage is CASAVA demultiplexing without index mismatches; PEAR merging at ≥50-base overlap and Q20; retention when ≥90% of bases exceed Q30; FASTX conversion; sed primer trimming; USEARCH zOTU processing; and BLASTn annotation. Exact versions, command lines, database snapshots, denoising/chimera settings and the control-exclusion rule have not been recovered as a complete executable workflow.

The deposit starts at retained counts, annotations and representatives. No FASTQs are included. Public raw-read accessions and their specimen linkage have not been verified. CP070274.1 is a third-party reference, not this study's raw-read accession.

## Software and reproduction

CSV, FASTA, Newick, Nexus, JSON and Markdown are open text formats. Source XLSX/XLSM requires compatible spreadsheet software; the CSVs require neither Excel nor macros. Worksheet exports contain cached stored values: formulas are not recalculated and macros are not executed.

The accompanying [Argiope_microbiome code repository](https://github.com/sven9r/Argiope_microbiome) documents the Makefile, packages and original project paths. R dependencies are recorded there in `renv.lock`; selected-run session records give more specific evidence where available. A current installed version is not automatically the version used for a historical figure. The code environment is not a complete FASTQ-processing environment.

Restore deposited inputs to their expected project locations before running analysis code:

| Deposit file | Code-project destination |
| --- | --- |
| `source_workbooks/Raw_data_COI_16S_2.xlsm` | `data/Raw_data_COI_16S_2.xlsm` |
| `source_workbooks/Argiope_Metadaten.xlsx` | `data/raw/Argiope_Metadaten.xlsx` |
| `source_workbooks/Tim_TabCOI_16S_1222.xlsx` | `data/raw/Tim_TabCOI_16S_1222.xlsx` |
| `third_party/host_phylogeny/` | `output/species_symbionts/host_phylogeny/` |
| `third_party/dusa_reference/` | `output/species_symbionts/reference_validation/` |

Original count/FASTA and saved-support paths are recorded in `file_manifest.csv`. Copying a file under a new name does not make every script runnable from the deposit root. Run commands from the **code-project root**, following its data-placement instructions.

| Make target | Purpose |
| --- | --- |
| `make restore` | Restore the R package environment recorded by the code repository. |
| `make saved-figures` | Retain source figures for the map and processing/collection-period/species-stack panels. |
| `make bacteria` | Bacterial profile states, focal log-ratio and dispersion outputs. |
| `make prey` | Prey composition and paired distance comparison. |
| `make history` | Historical DUSA-like tables and timeline selection. |
| `make variants` | Focal OTU variant composition. |
| `make composition` | Specimen-aligned bacterial/prey composition. |
| `make community-comparison` | Paired ordinations and clear-profile prey comparisons. |
| `make host-phylogeny` | Published host topology with local screening/composition summaries. |
| `make sequence-distance` | Five-fragment comparison and pairwise distances. |
| `make rarefaction` | Current read-depth/rarefaction display. |
| `make all` | Run the defined analysis/figure targets; consult the code README for dependencies and saved-input requirements. |

Rebuilding a scientific figure and copying a retained image are distinct operations. Stochastic calculations can vary unless inputs, parameters, random-number settings and package versions match. Use selected saved objects and session records when exact numerical correspondence matters. The figure map identifies where current-image sidecars or complete historical methods remain unavailable.

## Missing values, units and reuse cautions

- Empty CSV cells represent source blanks unless a file-specific definition says otherwise. Literal `NA` denotes missing/unavailable information. Existing tokens such as `NaN` can denote undefined fractions; the dictionary records tokens by column.
- Zero reads are observed zeros in an available count column. They differ from a missing feature column, excluded profile or missing metadata.
- Counts are integer reads. Fractions are dimensionless 0–1; explicitly named percentages use 0–100. Consult the dictionary rather than assuming that every generic `relative_abundance` field uses the same scale.
- Lengths and coordinates use bases. Stated fragment/gene coordinates are one-based and inclusive unless source provenance says otherwise. Newick branch units follow each tree's provenance.
- NMDS coordinates and Bray–Curtis dissimilarities are dimensionless. Ordination axes can rotate or reflect between runs without changing fitted distances.
- A zero or failing screen does not establish biological absence. A passing screen does not establish infection or exclude prey, environmental, laboratory or preservation sources.
- Geographic, host and temporal coverage is uneven. One contemporary collection does not represent the full species range.
- Provisional taxonomy is retained. The raw hit-label export lacks a complete accession/identity/coverage trail for independent validation of every label.
- Relative reads are compositional. These data do not measure absolute bacterial load, prey biomass, symbiont function, transmission or reproductive effects.
- Historical reads do not date acquisition or demonstrate uninterrupted infection. Recycled preservation ethanol is reported context, not documented here as the cause of a particular observation.

## Third-party material, citation and versioning

Cheng and Kuntner's host-tree dataset is [doi:10.5061/dryad.3c025](https://doi.org/10.5061/dryad.3c025), linked to [doi:10.1111/evo.12504](https://doi.org/10.1111/evo.12504). Its saved metadata record CC0. Reference genome CP070274.1 is associated in retained provenance with Busch et al. (2026), [doi:10.1093/femsle/fnag103](https://doi.org/10.1093/femsle/fnag103). Attribution and source terms are detailed in `THIRD_PARTY_SOURCES.md`.

The host display places 26 local species and leaves eight unplaced. Museum/contemporary *A. bruennichi* are separate sampling rows at one species tip. The plot does not estimate phylogenetic signal, ancestral infections or bacterial acquisition dates. Source branch lengths are retained in trees but suppressed in the topology display.

The current curated workbook's SHA-256 is `c9d1ca2992f83ab6e583ff00744e3f4d5f17217990312bee47079a707ca574f6`. The manifest identifies exact deposited bytes and excludes its own checksum to avoid circular hashing. Earlier workbooks and results remain dated provenance, not interchangeable versions of the current analysis.

Final creators, dataset citation, Dryad DOI, linked publication and project data licence must be entered at deposition. Source-specific CC0 does not license the new *Argiope* data. This documentation does not imply that a public release or submission has occurred.
