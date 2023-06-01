# Differential splicing detection

## Necessary packages/tools
- python
- R
- julia
- majiq/voila
- leafcutter (both leafcutter scripts & R package)
- whippet
- snakemake
- samtools
- regtools
- (R package) optparse
- (R package) data.table
- (R package) dplyr

### Leafcutter installation
Github installation: For leafcutter quantifications (junction counts for each splicing cluster), you will require the python scripts in `scripts` and `clustering` directories of the leafcutter github:
```
git clone https://github.com/davidaknowles/leafcutter
```
R package installation: For differential splicing analysis, the leafcutter R package is required:
```
devtools::install_github("davidaknowles/leafcutter/leafcutter")
```
Further detailed installation instructions can be found [here](https://davidaknowles.github.io/leafcutter/articles/Installation.html).

## Setting up

1. Upload the following 2 files in the `input` directory as follows:
`design.csv` : experimental design (change to your own, with the following columns)

| sample      | design    | genome | strand  | fq1             | fq2             |
|-------------|-----------|--------|---------|-----------------|-----------------|
| CTX_120 | control   | mm39   | reverse | CTX_120_1.fq.gz | CTX_120_2.fq.gz |
| CTX_125     | control   | mm39   | reverse | CTX_125_1.fq.gz | CTX_125_2.fq.gz |
| CTX_147     | control   | mm39   | reverse | CTX_147_1.fq.gz | CTX_147_2.fq.gz |
| CTX_148     | control   | mm39   | reverse | CTX_148_1.fq.gz | CTX_148_2.fq.gz |
| CTX_104 | treatment | mm39   | reverse | CTX_104_1.fq.gz | CTX_104_2.fq.gz |
| CTX_108     | treatment | mm39   | reverse | CTX_108_1.fq.gz | CTX_108_2.fq.gz |
| CTX_128     | treatment | mm39   | reverse | CTX_128_1.fq.gz | CTX_128_2.fq.gz |
| CTX_154     | treatment | mm39   | reverse | CTX_154_1.fq.gz | CTX_154_2.fq.gz |

`mouse.gff3` : gene annotation file (e.g. `gencode.vM29.primary_assembly.annotation.gff3`) You may obtain the annotation file from Gencode ([mouse](https://www.gencodegenes.org/mouse/), [human](https://www.gencodegenes.org/human/)).

2. Open the `prep.R` file and change the following parameters to suit your experimental design:

```
setwd("/mnt/cbis/home/yongshan/majiq_snakemake") # majiq snakemake directory
design <- read.csv("./input/design.csv") # experimental design csv file
bamdir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR" # directory containing your bam alignment files
experiment_name <- "WT_IgG2A_WT_O9_CTX" # set your own experiment name
gff3 <- "gencode.vM29.primary_assembly.annotation.gff3" # gff3 file name of species of interest

setwd("/mnt/cbis/home/yongshan/whippet_snakemake") # whippet snakemake directory
fasta_file_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa.gz" # location of genome fa.gz file
annotation_gtf_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz" # location of annotation gtf.gz file
experiment_name <- "WT_IgG2A_WT_O9_CTX" # set your own experiment name
design <- read.csv("./input/design.csv") # experimental design csv file
fq_dir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/trimmedFASTQ/" # directory with fastq files
bam_dir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR/" # directory with bam files
julia <- "/mnt/cbis/home/yongshan/julia-1.7.2/bin/julia" # location of julia command
whippet_bin <- "/mnt/cbis/home/yongshan/Whippet.jl/bin/" # directory of whippet scripts

setwd("/mnt/cbis/home/yongshan/leafcutter_snakemake") # leafcutter snakemake directory
experiment_name <- "WT_IgG2A_WT_O9_CTX" # set your own experiment name
design <- read.csv("./input/design.csv") # experimental design csv file
regtools_strand <- "1" # Note that: 0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR
gene_annotation <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz" # location of annotation gtf.gz file
bam_dir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR/" # directory with bam files
leafcutter_dir <- "/mnt/cbis/home/yongshan/leafcutter/" # directory of leafcutter installation
```

