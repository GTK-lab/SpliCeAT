import pandas as pd

sample_design = pd.read_table(config["SAMPLES"]).set_index("sample", drop=False)

SAMPLES = sample_design["sample"].tolist()

def get_bam_files(wildcards):
    return {sample_design.loc[wildcards.sample, "bam"]}