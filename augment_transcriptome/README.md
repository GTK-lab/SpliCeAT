# Generation of augmented transcriptome
**Use:** Augments a reference transcriptome with transcripts containing differential splicing events detected in the previous `ds_detection` step

**Input:** Differential splicing results from Majiq, Whippet & Leafcutter

**Output:** A kallisto index of the augmented transcriptome & transcript-to-gene (t2g) mapping for downstream quantification 

## Packages required
Click on the links for installation instructions. (In brackets are versions used when run by the author)
- [gffread](https://github.com/gpertea/gffread#installation) (0.12.7)
- [kallisto](https://pachterlab.github.io/kallisto/download) (0.46.1)
- [stringtie](https://ccb.jhu.edu/software/stringtie/#install) (2.2.1)

R packages:
- [yaml](https://www.rdocumentation.org/packages/yaml/versions/2.3.7) (2.3)
- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) (2.54.0)
- [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html) (1.5.0)
- [data.table](https://github.com/Rdatatable/data.table#installation) (1.14.6)
- [VennDetail](https://www.bioconductor.org/packages/release/bioc/html/VennDetail.html) (1.14.0)
- [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html) (1.1.0)
- [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) (1.50.2)

## Required input
The `ds_detection` pipeline should have been run and the following output files should be present: 
1. MAJIQ `LSV` file (e.g. `lsvs.txt`)
2. Whippet `.diff` file (e.g. `tdp43_nestin_ctx_e14.diff`)
3. Leafcutter `cluster_significance.txt` and `effect_sizes.txt` files (e.g. `tdp43_nestin_ctx_e14_cluster_significance.txt` and `tdp43_nestin_ctx_e14_effect_sizes.txt`)

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

2. Open `config/config.yaml` and change the following paths:
```
BASE_PATH: /mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/
SAMPLES: /mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/config/samples.tsv

###################
# STRINGTIE ASSEMBLY
###################
# If you have added the executable to your PATH, simply put "stringtie" for STRINGTIE_COMMAND:
STRINGTIE_COMMAND: /mnt/cbis/home/yongshan/stringtie-2.2.1.Linux_x86_64/stringtie
GTF: /mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf
STRINGTIE_THREADS: 4

# minimum junction coverage for Stringtie assembly:
JUNCS_CUTOFF: 20 

# --rf : assume stranded library fr-firststrand; --fr : assume stranded library fr-secondstrand:
STRAND: rf

# alignments directory should contain all sample BAM & BAI files:
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
# If you have added the executable to your PATH, simply put "gffread" for gffread_path:
gffread_path: /mnt/cbis/home/yongshan/gffread/gffread
genome_fasta: /mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa
transcripts_fasta: /mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.transcripts.fa
```

3. Open up `prep_for_AT.R` & `collapse_transcripts.R` in the `scripts` folder & the `Snakefile` in `workflow` and change the first line to point to your config file location:

`prep_for_AT.R`/`collapse_transcripts.R`:
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
Once the above files are configured correctly, execute a Snakemake dry run with
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
|   ├── merged_stringtie_assembly_novel_exon_filtered.fa
|   ├── t2g_augment_uncollapsed.csv
│   └── t2g_augment_collapsed.csv
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
The transcriptome index is found at `kallisto_index_augmented_transcriptome` and can be used for subsequent kallisto quantification and sleuth differential expression analysis. The t2g (both uncollapsed `t2g_augment_uncollapsed.csv` and collapsed `t2g_augment_collapsed.csv` versions) mappings can be used for subsequent differential expression analysis. 
