# Majiq snakemake

## Useful tools
- [Majiq command builder](https://biociphers.bitbucket.io/majiq-docs-academic/commandbuilder.html)

## Necessary packages/tools
- snakemake
- majiq/voila
- dplyr

## Setting up
1. Open the `prep.R` file and change the following parameters to suit your experimental design:
```
setwd("/mnt/cbis/home/yongshan/majiq_snakemake") # majiq snakemake directory
design <- read.csv("./input/design.csv") # experimental design csv file, do not need to change if you upload it under input
bamdir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR" # directory containing your bam alignment files
experiment_name <- "WT_IgG2A_WT_O9_CTX" # set your own experiment name
gff3 <- "gencode.vM29.primary_assembly.annotation.gff3" # gff3 file name of species of interest
base_path <- getwd()
```

2. Open the `Snakefile` in `workflow` and change the first line to point to your config file location:
```
 configfile: "<your_snakemake_dir>/config/config.yaml"
 ...
 ```

3. Upload the following 2 files in the `input` directory as follows:

`design.csv` : experimental design (change to your own, with the following columns)

| sample | design    | genome | strand  |
| ------ | ------    | ------ | ------  |
| F2014	 | control	 | mm39	  | reverse |
| F2011	 | control	 | mm39	  | reverse |
| F2003	 | control	 | mm39	  | reverse |
| F2001	 | treatment | mm39	  | reverse |
| F2010	 | treatment | mm39	  | reverse |
| F2002	 | treatment | mm39	  | reverse |

`mouse.gff3` : gene annotation file (e.g. `gencode.vM29.primary_assembly.annotation.gff3`) You may obtain the annotation file from Gencode ([mouse](https://www.gencodegenes.org/mouse/), [human](https://www.gencodegenes.org/human/)).

4. Run `prep.R` on command line with
```
Rscript prep.R
```
in order to populate the directories with the necessary helper files. The following files should be successfully created:
- `conf/<experiment_name>_conf.txt`
- `config/config.yaml`
- `config/confs.tsv`
- `config/delta_psi_samples.tsv`
- `config/experiment_sample_names.tsv`


## Running Snakemake

Once the above finishes running successfully and the necessary helper files are created, execute a Snakemake dry run with
```
snakemake -np
```
to check the parameters of the run. Once ready to run, execute
```
snakemake --cores 24
```

## Output

In the `results` directory, you should see the following output files:
`results/majiq_build`:
- `<sample>.majiq`
- `<sample>.sj`
- `majiq.log`
- `splicegraph.sql`

`results/majiq_delta_psi`:
- `<group_comparison>.deltapsi.tsv`
- `<group_comparison>.deltapsi.voila`
- `deltapsi_majiq.log`

In order to view the interactive results on Voila web server, execute the following on command line:
```
cd results/majiq_delta_psi/<experiment_name>
voila view <group_comparison>.deltapsi.voila results/majiq_build/<experiment_name>/splicegraph.sql -p 5000 --host 0.0.0.0
```
View the results on http://<remote_server>:5000/ (e.g. `http://atlas.cbis.nus.edu.sg:5000/`) and they should be displayed as follows:
![image](https://github.com/ys-lim/SpliCeAT/assets/68455070/83af4e97-5fd7-440b-9f33-28bd2fecc5a2)
Click on `Download LSVs` to get the list of local splicing variations. 
