import pandas as pd

samples = pd.read_table(config["samples_tsv"]).set_index("experiment", drop=False)
fastq = pd.read_table(config["fastq_tsv"]).set_index("sample", drop=False)
EXPERIMENTS = samples["experiment"].tolist()
SAMPLES = fastq["sample"].tolist()
delta = pd.read_table(config["delta_tsv"]).set_index("experiment", drop=False)
delta_input_tsv = pd.read_table(config["delta_input_tsv"]).set_index("experiment", drop=False)

def samtools_input(wildcards):
    string = samples.loc[wildcards.experiment, "cko"]
    string_list = string.split()
    new_string_list = []
    for i in string_list:
        new_i = samples.loc[wildcards.experiment, "dir"]+i
        new_string_list.append(new_i)
    return new_string_list
	
def get_fastq(wildcards):
    fq1 = fastq.loc[wildcards.sample, "dir"]+fastq.loc[wildcards.sample, "fq1"]
    fq2 = fastq.loc[wildcards.sample, "dir"]+fastq.loc[wildcards.sample, "fq2"]
    return [fq1,fq2]
	
def experiment_sample(wildcards):
    output = config["base_path"]+"/results/quantify/"+fastq.loc[wildcards.experiment, "experiment"]+"/"+fastq.loc[wildcards.experiment, "sample"]+"/"+fastq.loc[wildcards.experiment, "sample"]
    return output
	
def get_index(wildcards):
    return fastq.loc[wildcards.sample, "index_dir"]
	
def grp1(wildcards):
    return delta.loc[wildcards.experiment, "grp1"]
	
def grp2(wildcards):
    return delta.loc[wildcards.experiment, "grp2"]
	
def delta_input(wildcards):
    return delta_input_tsv.loc[wildcards.experiment, "psi_gz"]