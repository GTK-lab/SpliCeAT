# Generation of augmented transcriptome
This pipeline augments a reference transcriptome with transcripts containing differential splicing events detected in the previous `ds_detection` step. The output of this pipeline is a kallisto index of the augmented transcriptome.

## Required input
The `ds_detection` pipeline should have been run and the following files should be present: 
- 
These files are used as input for generating the augmented transcriptome.
