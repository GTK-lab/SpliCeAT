# Generation of augmented transcriptome
This pipeline augments a reference transcriptome with transcripts containing differential splicing events detected in the previous `ds_detection` step. The output of this pipeline is a `kallisto` index of the augmented transcriptome.

## Required input
The `ds_detection` pipeline should have been run and the following output files should be present: 
1. MAJIQ `LSV` file 
2. Whippet `.diff` file
3. Leafcutter `cluster_significance.txt` and `effect_sizes.txt` files

These files contain high-confidence differential splicing events which will subsequently be included in the augmented transcriptome.
