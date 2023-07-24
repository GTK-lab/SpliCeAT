# Differential expression analysis

**Use:** Performs differential transcript expression analysis using augmented transcriptome (generated in previous [step](https://github.com/ys-lim/SpliCeAT/tree/main/augment_transcriptome))

**Input:** kallisto index of the augmented transcriptome & sample FASTQ files

**Output:** List of differentially expressed transcripts (in both uncollapsed and collapsed t2g versions)

## Packages required
Click on the links for installation instructions. (In brackets are versions used when run by the author)
- [kallisto](https://pachterlab.github.io/kallisto/download) (0.46.1)
- [sleuth](http://pachterlab.github.io/sleuth/download) (0.30.1)

R packages:
- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) (2.54.0)
- [lgr](https://cran.r-project.org/web/packages/lgr/index.html) (0.4.4)

## Required input
The `ds_detection` and `augment_transcriptome` steps should have been run and the following output files (located in `SpliCeAT/augment_transcriptome/`) will be used in this step: 
```
results
└── augmented_transcriptome
    ├── kallisto_index_augmented_transcriptome
    ├── t2g_augment_collapsed.csv
    └── t2g_augment_uncollapsed.csv
```
The subsequent analysis will be carried out in 2 parts: (1) Quantification of transcripts against the augmented transcriptome, and (2) Differential transcript expression analysis. 

## Setting up
