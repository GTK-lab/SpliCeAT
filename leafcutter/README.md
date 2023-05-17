# Leafcutter snakemake

## Necessary packages/tools
- leafcutter (both leafcutter scripts & R package)
- optparse
- dplyr
- data.table
- regtools

### Leafcutter installation

For quantification:
```
git clone https://github.com/davidaknowles/leafcutter
```
For leafcutter quantifications (junction counts for each splicing cluster), you will require the python scripts in `scripts` and `clustering` directories.

For differential splicing analysis, the leafcutter R package is required:
```
devtools::install_github("davidaknowles/leafcutter/leafcutter")
```
Detailed installation instructions are [here](https://davidaknowles.github.io/leafcutter/articles/Installation.html).

## Setting up

1. Upload the `design.csv` experimental design file into the `input` directory as follows:

| sample | design  |
| ------ | ------  |
| F2014	 | control |
| F2011	 | control |
| F2003	 | control |
| F2001	 | treatment |
| F2010	 | treatment |
| F2002	 | treatment |

2. Open the `prep.R` file and change the following parameters to suit your experimental design:
```
setwd("/mnt/cbis/home/yongshan/leafcutter_snakemake") # leafcutter snakemake directory
base_path <- getwd()
experiment_name <- "WT_IgG2A_WT_O9_CTX" # set your own experiment name
design <- read.csv("./input/design.csv") # experimental design csv file
regtools_strand <- "1" # Note that: 0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR
gene_annotation <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz" # location of annotation gtf.gz file
bam_dir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR/" # directory with bam files
leafcutter_dir <- "/mnt/cbis/home/yongshan/leafcutter/" # directory of leafcutter
```

## Setting up
1. ./gtf_to_exons.R /mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz ./test.txt.gz

## Other useful tools
- [Leafcutter vignette](https://davidaknowles.github.io/leafcutter/articles/Usage.html)

