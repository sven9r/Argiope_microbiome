# Reference data

## Spider phylogeny

Figure 2 uses the Bayesian host topology `Bestree_BI10M` from Cheng and Kuntner (2014), *Phylogeny suggests nondirectional and isometric evolution of sexual size dimorphism in argiopine spiders* ([article](https://doi.org/10.1111/evo.12504); [Dryad dataset](https://doi.org/10.5061/dryad.3c025)). The deposited data are identified as CC0 in the retained source metadata.

The original Dryad `Appendix S2` file has MD5 `c608bba1562615ceb4ea07bb45b4b5e8`. The retained derivative `output/species_symbionts/host_phylogeny/host_species_tree_for_plot.nwk` contains 26 matched species. `species_tip_mapping.csv` records all 34 local species, including the eight without a matching tip. The displayed tree shows topology; it does not show divergence dates or bacterial transmission. Original host loci were cox1, cox2, H3 and 28S.

## Deposited bacterial reference

Supplementary Figure 8 uses the annotated 16S rRNA gene from [GenBank CP070274.1](https://www.ncbi.nlm.nih.gov/nuccore/CP070274.1). The saved GenBank record was retrieved on 7 October 2026. Its SHA256 is `821a56383a27eb78a319460a56a195fe8e69525cc9f404f26a5d738473046fe6`.

The feature occupies genomic positions 553873–555402 on the minus strand (one-based inclusive coordinates). `CP070274_16S.fasta` contains the oriented 1,530-base gene. The renderer checks that extraction against the saved record. Local OTU2 matches positions 1060–1160 of the oriented gene at all 101 measured bases. This short-fragment match does not establish strain identity or infection status.

These source attributions concern the reference records. No project-wide code or data licence has been assigned to this private preparation repository.
