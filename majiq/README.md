setwd("/mnt/cbis/home/yongshan/majiq_snakemake")
design <- read.csv("./input/design.csv")
bamdir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR"
experiment_name <- "WT_IgG2A_WT_O9_CTX"
gff3 <- "gencode.vM29.primary_assembly.annotation.gff3"
base_path <- getwd()
