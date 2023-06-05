import pandas as pd

experimental_design_majiq = pd.read_table(config["confs"]).set_index("experiment", drop=False)
samples_majiq = pd.read_table(config["delta_psi_samples"]).set_index("experiment", drop=False)
experiment_samples_majiq = pd.read_table(config["experiment_sample_names"]).set_index("experiment", drop=False)
EXPERIMENTS_majiq = experimental_design["experiment"].tolist()
SAMPLE_NAMES_majiq = experiment_samples["experiment"].tolist()

def get_conf_file(wildcards):
    return{
        experimental_design_majiq.loc[wildcards.experiment, "conf_path"]
    }
	
def delta_psi_grp1(wildcards):
    string = samples_majiq.loc[wildcards.experiment, "grp1"]
    string_list = string.split()
    new_string_list = []
    for i in string_list:
        new_i = samples_majiq.loc[wildcards.experiment, "dir"]+i+".majiq"
        new_string_list.append(new_i)
    return new_string_list

def delta_psi_grp2(wildcards):
  string = samples_majiq.loc[wildcards.experiment, "grp2"]
  string_list = string.split()
  new_string_list = []
  for i in string_list:
      new_i = samples_majiq.loc[wildcards.experiment, "dir"]+i+".majiq"
      new_string_list.append(new_i)
  return new_string_list
  
  
samples_leafcutter = pd.read_table(config["samples_tsv"]).set_index("sample", drop=False)
juncs_leafcutter = pd.read_table(config["juncs_file"]).set_index("experiment", drop=False)
output_junc_leafcutter = pd.read_table(config["output_juncs"]).set_index("experiment", drop=False)
SAMPLES_leafcutter = samples["sample"].tolist()
EXPERIMENTS_leafcutter = juncs["experiment"].tolist()

def get_bam(wildcards):
    output = config["bam_dir"]+samples_leafcutter.loc[wildcards.sample, "bam"]
    return output

def get_junc(wildcards):
    return config["BASE_PATH"]+"/config/"+juncs_leafcutter.loc[wildcards.experiment,"juncs_file"]
	
def get_output_junc(wildcards):
    return config["BASE_PATH"]+"/results/"+output_junc_leafcutter.loc[wildcards.experiment,"junction_files"]


samples_whippet = pd.read_table(config["samples_tsv"]).set_index("experiment", drop=False)
fastq_whippet = pd.read_table(config["fastq_tsv"]).set_index("sample", drop=False)
EXPERIMENTS_whippet = samples["experiment"].tolist()
SAMPLES_whippet = fastq["sample"].tolist()
delta_whippet = pd.read_table(config["delta_tsv"]).set_index("experiment", drop=False)
delta_input_tsv_whippet = pd.read_table(config["delta_input_tsv"]).set_index("experiment", drop=False)

def samtools_input(wildcards):
    string = samples_whippet.loc[wildcards.experiment, "cko"]
    string_list = string.split()
    new_string_list = []
    for i in string_list:
        new_i = samples_whippet.loc[wildcards.experiment, "dir"]+i
        new_string_list.append(new_i)
    return new_string_list
	
def get_fastq(wildcards):
    fq1 = fastq_whippet.loc[wildcards.sample, "dir"]+fastq_whippet.loc[wildcards.sample, "fq1"]
    fq2 = fastq_whippet.loc[wildcards.sample, "dir"]+fastq_whippet.loc[wildcards.sample, "fq2"]
    return [fq1,fq2]
	
def experiment_sample(wildcards):
    output = config["BASE_PATH"]+"/results/quantify/"+fastq_whippet.loc[wildcards.experiment, "experiment"]+"/"+fastq_whippet.loc[wildcards.experiment, "sample"]+"/"+fastq_whippet.loc[wildcards.experiment, "sample"]
    return output
	
def get_index(wildcards):
    return fastq_whippet.loc[wildcards.sample, "index_dir"]
	
def grp1(wildcards):
    return delta_whippet.loc[wildcards.experiment, "grp1"]
	
def grp2(wildcards):
    return delta_whippet.loc[wildcards.experiment, "grp2"]
	
def delta_input(wildcards):
    return delta_input_tsv_whippet.loc[wildcards.experiment, "psi_gz"]
