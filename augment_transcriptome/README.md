# Generation of augmented transcriptome
This pipeline augments a reference transcriptome with transcripts containing differential splicing events detected in the previous `ds_detection` step. The output of this pipeline is a `kallisto` index of the augmented transcriptome.

## Necessary packages/tools

## Required input
The `ds_detection` pipeline should have been run and the following output files should be present: 
1. MAJIQ `LSV` file 
2. Whippet `.diff` file
3. Leafcutter `cluster_significance.txt` and `effect_sizes.txt` files

These files contain high-confidence differential splicing events which will subsequently be included in the augmented transcriptome.

## Setting up
1. Open the `samples.tsv` in the `config` folder and change to your sample names accordingly:
```
sample
CTX_120
CTX_125
CTX_147
CTX_148
CTX_104
CTX_108
CTX_128
CTX_154
```

2. Open up `config.yaml` and change the following paths:
```
BASE_PATH: /mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/
SAMPLES: /mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/config/samples.tsv

###################
# STRINGTIE ASSEMBLY
###################
STRINGTIE_COMMAND: /mnt/cbis/home/yongshan/stringtie-2.2.1.Linux_x86_64/stringtie
GTF: /mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf
# minimum junction coverage for Stringtie assembly
JUNCS_CUTOFF: 20 
STRINGTIE_THREADS: 4
# --rf : assume stranded library fr-firststrand; --fr : assume stranded library fr-secondstrand
STRAND: rf
ALIGNMENTS_DIR: /mnt/gtklab01/linglab/tdp43/STAR/tdp43_nestin_ctx_e14/

###################
# MASTERLISTS PREP
###################
mouse_ensembl_version: 106
majiq_lsv_file_path: /mnt/cbis/home/yongshan/SpliCeAT/ds_detection/results/majiq_delta_psi/tdp43_nestin_ctx_e14/lsvs.txt
leafcutter_cluster_sig_file_path: /mnt/cbis/home/yongshan/SpliCeAT/ds_detection/results/tdp43_nestin_ctx_e14_cluster_significance.txt
leafcutter_effect_size_file_path: /mnt/cbis/home/yongshan/SpliCeAT/ds_detection/results/tdp43_nestin_ctx_e14_effect_sizes.txt
whippet_diff_file_path: /mnt/cbis/home/yongshan/SpliCeAT/ds_detection/results/delta_psi/tdp43_nestin_ctx_e14.diff

###################
# AUG TRANSCRIPTOME
###################
gffread_path: /mnt/cbis/home/yongshan/gffread/gffread
genome_fasta: /mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa
transcripts_fasta: /mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.transcripts.fa
```

3. Open up `prep_for_AT.R` in the `scripts` folder & the `Snakefile` in `workflow` and change the first line to point to your config file location:

`prep_for_AT.R`:
```
#### CHANGE THIS ####
config_file_path <- "/mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/config/config.yaml"
####################
...
```
`Snakefile`:
```
configfile: "/mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/config/config.yaml"
...
```

## Running Snakemake
Once the above finishes running successfully and the necessary helper files are created, execute a Snakemake dry run with
```
snakemake -np
```
to check the parameters of the run. Once ready to run, execute
```
snakemake --cores 24
```

## Output files
