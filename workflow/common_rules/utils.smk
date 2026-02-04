from snakemake.utils import validate
import pandas as pd
import os
import yaml
import shutil
from pathlib import Path

PROJECT_ROOT = config["pipeline_dir"]
results_base = config.get("results_dir") or PROJECT_ROOT
RESULTS_DIR = os.path.join(results_base, "results")

# result directories
R0_REF = os.path.join(RESULTS_DIR, "r0_get_ref")
R1_DS = os.path.join(RESULTS_DIR, "r1_ds_detection")
R2_AT = os.path.join(RESULTS_DIR, "r2_augment_transcriptome")
R3_DE = os.path.join(RESULTS_DIR, "r3_de_analysis")

# for reference files
RF_DIR = os.path.join(R0_REF, config['ref']['species'])
STR_DIR = os.path.join(R0_REF, "STAR")

# for ds_detection
LC_DIR = os.path.join(R1_DS, "leafcutter")
RT_DIR = os.path.join(R1_DS, "regtools")
MJ_DIR = os.path.join(R1_DS, "majiq")
WP_DIR = os.path.join(R1_DS, "whippet")

# for de_analysis
KQ_DIR = os.path.join(R3_DE, "kallisto_quant_out")
SL_DIR = os.path.join(R3_DE, "sleuth")

# config resolve pathnames
original_design_path = os.path.join(PROJECT_ROOT,'config',config["samples"])
original_df=pd.read_csv(original_design_path, sep="\t", dtype=str, comment="#")

star_active = config.get('STAR', {}).get('activate', False)
samples_full_path = os.path.join(RESULTS_DIR, "samples.tsv")

if os.path.exists(samples_full_path):
	sample_file_df = pd.read_csv(samples_full_path, sep="\t", dtype=str, comment="#")
else:
	sample_file_df = original_df.copy()
	if star_active and 'bam_file' not in sample_file_df.columns:
		sample_file_df['bam_file'] = [
			os.path.join(STR_DIR, "aligned_BAM", f"{s}_Aligned.sortedByCoord.out.bam")
			for s in sample_file_df['sample_name']
		]

	if not star_active and 'bam_file' not in sample_file_df.columns:
		raise ValueError(
		"\n[Metadata Error]: Inonsistent input files- STAR is disabled ('activate: false') in config.yaml, but design.tsv is missing the required 'bam_file' column. Please provide paths or enable STAR."
		)

	if star_active and 'bam_file' in sample_file_df.columns:
		raise ValueError(
		"\n[Metadata Error]: Inonsistent input files- STAR is enabled ('activate: true'), but the 'bam_file' column already exists in design.tsv. To use your own BAMs please disable STAR in the config. To use STAR, remove the 'bam_file' column."
	)

	os.makedirs(RESULTS_DIR, exist_ok=True)
	sample_file_df.to_csv(samples_full_path, sep="\t", index=False)

annot = sample_file_df.set_index("sample_name", drop=False)

annot['bam_dirs']= [str(f.parent) for f in [ Path(f) for f in annot['bam_file'] ]]
annot['bam_stem'] = [str(f.stem) for f in [ Path(f) for f in annot['bam_file'] ]]  # .stem removes the extension

SAMPLES = annot["sample_name"].tolist() # list of all sample names
GROUPS = config["experiment"]["groups"]
# list of comparison groups

# convert annot into dictionary for parametrization of rules, by deduplicating by sample_name (should only differ by bam_file)
samples = annot.to_dict(orient="index")

bam_files = list(annot['bam_file'])

def get_sample_name(wildcards):
	return annot.loc[wildcards.sample, "sample_name"]

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

ds_det_output = [
	os.path.join(MJ_DIR,f"majiq_expanded_{'-'.join(GROUPS)}.deltapsi.tsv"),
	os.path.join(WP_DIR, "whippet_delta_psi.diff"),
	os.path.join(LC_DIR, "leafcutter_lsvs.tsv")
]

# for augment_transcriptome
AT_DIR = os.path.join(R2_AT, "augmented_transcriptome")
MG_DIR = os.path.join(R2_AT, "merged_assembly")
ST_DIR = os.path.join(R2_AT, "stringtie_assemblies")
ML_DIR = os.path.join(R2_AT, "masterlists")


