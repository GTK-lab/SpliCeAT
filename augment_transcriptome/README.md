# Generation of augmented transcriptome
This pipeline augments a reference transcriptome with transcripts containing differential splicing events detected in the previous `ds_detection` step. The output of this pipeline is a `kallisto` index of the augmented transcriptome.

## Necessary packages/tools
- [gffread](https://github.com/gpertea/gffread#installation) (0.12.7)
- [kallisto](https://pachterlab.github.io/kallisto/download) (0.46.1)
- [stringtie](https://ccb.jhu.edu/software/stringtie/#install) (2.2.1)
- (R package) yaml (2.3)
- (R package) biomaRt (2.54.0)
- (R package) stringr (1.5.0)
- (R package) data.table (1.14.6)
- (R package) VennDetail (1.14.0)
- (R package) dplyr (1.1.0)
- (R package) GenomicRanges (1.50.2)

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
# If you have added the executable to your PATH, simply put "stringtie"
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
# If you have exported to path, simply put "gffread"
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
You should see the following files in the `results` folder:
```
results
├── augmented_transcriptome
│   ├── augmented_transcripts.fa
│   ├── kallisto_index_augmented_transcriptome
│   └── merged_stringtie_assembly_novel_exon_filtered.fa
├── merged_assembly
│   ├── merged_stringtie_assembly.gtf
│   ├── merged_stringtie_assembly_novel_exon_filtered.gtf
│   └── merged_stringtie_assembly_novel_exon_filtered_with_reference.gtf
└── stringtie_assemblies
    ├── <sample1>_ref_guided_assembly.gtf
    ├── <sample2>_ref_guided_assembly.gtf
    ├── <sample3>_ref_guided_assembly.gtf
    ├── <sample4>_ref_guided_assembly.gtf
    ├── <sample5>_ref_guided_assembly.gtf
    ├── <sample6>_ref_guided_assembly.gtf
    ├── <sample7>_ref_guided_assembly.gtf
    └── <sample8>_ref_guided_assembly.gtf
```
The transcriptome index is found at `kallisto_index_augmented_transcriptome` and will be used for subsequent kallisto quantification and sleuth differential expression analysis.
