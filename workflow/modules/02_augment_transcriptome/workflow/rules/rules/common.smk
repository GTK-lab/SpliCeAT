import pandas as pd

sample_file = pd.read_csv(samples_full_path, sep="\t", dtype=str, comment="#")
sample_design = sample_file.set_index("sample_name", drop=False)

SAMPLES = sample_design["sample_name"].tolist()

def get_bam(wildcards):
    return [sample_design.loc[wildcards.sample, "bam_file"]]