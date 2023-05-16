import pandas as pd

experimental_design = pd.read_table(config["confs"]).set_index("experiment", drop=False)
samples = pd.read_table(config["delta_psi_samples"]).set_index("experiment", drop=False)
experiment_samples = pd.read_table(config["experiment_sample_names"]).set_index("experiment", drop=False)

EXPERIMENTS = experimental_design["experiment"].tolist()
SAMPLE_NAMES = experiment_samples["experiment"].tolist()

def get_conf_file(wildcards):
    return{
        experimental_design.loc[wildcards.experiment, "conf_path"]
    }
	
def delta_psi_grp1(wildcards):
    string = samples.loc[wildcards.experiment, "grp1"]
    string_list = string.split()
    new_string_list = []
    for i in string_list:
        new_i = samples.loc[wildcards.experiment, "dir"]+i+".majiq"
        new_string_list.append(new_i)
    return new_string_list

def delta_psi_grp2(wildcards):
    string = samples.loc[wildcards.experiment, "grp2"]
    string_list = string.split()
    new_string_list = []
    for i in string_list:
        new_i = samples.loc[wildcards.experiment, "dir"]+i+".majiq"
        new_string_list.append(new_i)
    return new_string_list