from snakemake.utils import validate
import pandas as pd
import os
import yaml
from pathlib import Path

PROJECT_ROOT = config["pipeline_dir"]
RESULTS_DIR = config.get("results_dir", os.path.join(PROJECT_ROOT, "results"))

# config resolve pathnames
sample_file_name=config["samples"]
samples_full_path = os.path.join(PROJECT_ROOT,'config',sample_file_name)

sample_file_df = pd.read_csv(samples_full_path, sep="\t", dtype=str, comment="#")

annot = sample_file_df.set_index("sample_name", drop=False)

annot['bam_dirs']= [str(f.parent) for f in [ Path(f) for f in annot['bam_file'] ]]
annot['bam_stem'] = [str(f.stem) for f in [ Path(f) for f in annot['bam_file'] ]]  # .stem removes the extension

SAMPLES = annot["sample_name"].tolist() # list of all sample names
GROUPS = annot["group"].unique().tolist() if "group" in annot.columns else [] # list of comparison groups

# convert annot into dictionary for parametrization of rules, by deduplicating by sample_name (should only differ by bam_file)
samples = annot.to_dict(orient="index")

bam_files = list(annot['bam_file'])

# fasta file name lists
def get_fq1(wildcards):
	return annot.loc[wildcards.sample, "fq1"]

def get_fq2(wildcards):
	return annot.loc[wildcards.sample, "fq2"]

def get_fastq_files(wildcards):
	return {
		"read1": annot.loc[wildcards.sample, "fq1"],
		"read2": annot.loc[wildcards.sample, "fq2"]
	}

def get_bam(wildcards):
	return annot.loc[wildcards.sample, "bam_file"]

def get_bam_stem(wildcards):
	return annot.loc[wildcards.sample, "bam_stem"]

def get_group(wildcards):
	return annot.loc[wildcards.sample, "group"]

def sample_for_group(group):
	return list(annot[annot['group'] == group].index)

def bamfiles_for_group(group):
	return list(annot[annot['group'] == group]['bam_file'])

def bamfiles_for_groups(groups):
	return [bamfiles_for_group(group) for group in groups]

def bamstems_for_group(group):
	return list(annot[annot['group'] == group]['bam_stem'])

def bamstems_for_groups(groups):
	return [bamstems_for_group(group) for group in groups]