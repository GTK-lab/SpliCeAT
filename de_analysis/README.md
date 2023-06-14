# Differential expression analysis

## Necessary packages/tools
- kallisto (0.46.1)
- (R package) sleuth (0.30.1)

## Required input
The `ds_detection` and `augment_transcriptome` pipelines should have been run and the following output files should be present: 
- kallisto

This transcriptome index contains both reference transcripts and novel transcripts with differential splicing events. We will now test for differential expression at the gene and transcript level.
