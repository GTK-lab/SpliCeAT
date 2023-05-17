import pandas as pd

samples = pd.read_table(config["samples_tsv"]).set_index("sample", drop=False)
SAMPLES = samples["sample"].tolist()
juncs = pd.read_table(config["juncs_file"]).set_index("experiment", drop=False)
EXPERIMENTS = juncs["experiment"].tolist()
output_junc = pd.read_table(config["output_juncs"]).set_index("experiment", drop=False)

def get_bam(wildcards):
    output = config["bam_dir"]+samples.loc[wildcards.sample, "bam"]
    return output

def get_junc(wildcards):
    return config["BASE_PATH"]+"config/"+juncs.loc[wildcards.experiment,"juncs_file"]
	
def get_output_junc(wildcards):
    return config["BASE_PATH"]+"results/"+output_junc.loc[wildcards.experiment,"junction_files"]