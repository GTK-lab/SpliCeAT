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

1. Upload your `design.csv` in the `input` directory as follows:

`design.csv` : experimental design (change to your own, with the following compulsory columns)

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

2. Download the following gene annotations and genome files and place them into a location of your choice:

- `GFF3` gene annotation file  (e.g. `gencode.vM29.primary_assembly.annotation.gff3`).  Unzip the file using `gunzip` tool.
- `GTF` gene annotation file (e.g. `gencode.vM29.primary_assembly.annotation.gtf.gz`).
- `FASTA` genome file (e.g. `GRCm39.primary_assembly.genome.fa.gz`).

You may obtain the annotation files from Gencode ([mouse](https://www.gencodegenes.org/mouse/), [human](https://www.gencodegenes.org/human/)).

3. Open the `prep.R` file and change the following parameters to suit your experimental design:

```
# differential splicing detection directory
setwd("/mnt/cbis/home/yongshan/SpliCeAT/ds_detection") 

# experimental design csv file - dont need to change if design.csv is in input directory
design <- read.csv("./input/design.csv") 

# directory containing your bam alignment files
bam_dir <- "/mnt/gtklab01/linglab/tdp43/STAR/tdp43_nestin_ctx_e14/" 

 # directory containing your fastq files
fq_dir <- "/mnt/gtklab01/linglab/tdp43/fastq/"

# set your own experiment name
experiment_name <- "tdp43_nestin_ctx_e14" 

# path of annotation gff3 file
gff3_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gff3" 

# path of annotation gtf.gz file
annotation_gtf_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz" 

# path of genome fa.gz file
fasta_file_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa.gz" 

# path of julia command
julia <- "/mnt/cbis/home/yongshan/julia-1.7.2/bin/julia" 

# directory of whippet scripts
whippet_bin <- "/mnt/cbis/home/yongshan/Whippet.jl/bin/" 

# Note that: 0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR
regtools_strand <- "1" 

# directory of leafcutter installation
leafcutter_dir <- "/mnt/cbis/home/yongshan/leafcutter/" 
```

4. Open the `Snakefile` in `workflow` and change the first line to point to your config file location:
```
 configfile: "<your_ds_detection_snakemake_dir>/config/config.yaml"
 ...
 ```

5. Run `prep.R` on command line with
```
Rscript prep.R
```
in order to populate the directories with the necessary helper files. The following 14 files should be successfully created:
- `config/config.yaml`

### Majiq helper files
- `config/<experiment_name>_conf.txt`
- `config/confs.tsv`
- `config/delta_psi_samples.tsv`
- `config/experiment_sample_names.tsv`

### Whippet helper files
- `config/delta.tsv`
- `config/delta_input.tsv`
- `config/fastq.tsv`
- `config/samples_whippet.tsv`

### Leafcutter helper files
- `config/juncs_file.tsv`
- `config/output_junc.tsv`
- `config/samples_leafcutter.tsv`
- `config/<experiment_name>_groups_file.txt`
- `config/<experiment_name>_groups_junc.txt`
