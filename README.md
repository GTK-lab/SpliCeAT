<p align="center">
  <img src="images/Logo.png">
</p>

# SpliCeAT: Integrated pipeline for detection and quantification of aberrant transcripts with novel splicing events

This repository contains the following steps of the pipeline [^1]:
1. Majiq
2. Whippet
3. Leafcutter
4. Augmented transcriptome analysis

Start here: 
```
git clone https://github.com/ys-lim/SpliCeAT.git
```

[^1]: All RNA-seq alignments are performed using STAR. The pipelines expect BAM files to be labelled as `sample_Aligned.sortedByCoord.out.bam`. Nevertheless, modifications can be made (at the user's discretion) in the Snakemake `rules` to account for alignments generated by other tools (e.g. HISAT2). 
