# Majiq snakemake

## Setting up
1. In the `prep.R` file, change the following parameters to suit your experimental design:
`setwd("/mnt/cbis/home/yongshan/majiq_snakemake") # working directory containing majiq snakemake
design <- read.csv("./input/design.csv") # experimental design csv file, do not need to change if you upload it under input
bamdir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR" # directory containing your bam alignment files
experiment_name <- "WT_IgG2A_WT_O9_CTX" # set your own experiment name
gff3 <- "gencode.vM29.primary_assembly.annotation.gff3" # gff3 file name of species of interest
base_path <- getwd()`
2. Upload your experimental design csv file in `input` directory as follows:
| sample  | design | genome  | strand |
| ------- | ------ | ------ | ------- |
| F2014	| control	| mm39	| reverse |
| F2011	| control	| mm39	| reverse |
| F2003	| control	| mm39	| reverse |
| F2001	| treatment	| mm39	| reverse |
| F2010	| treatment	| mm39	| reverse |
| F2002	| treatment	| mm39	| reverse |

