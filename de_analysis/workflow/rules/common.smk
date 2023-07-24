import pandas as pd

experimental_design = pd.read_table(config["BASE_PATH"]+"/de_analysis/config/samples.tsv").set_index("sample_name", drop=False)

SAMPLES = experimental_design["sample_name"].tolist()

def get_fastq_files(wildcards):
    return{
        "read1": experimental_design.loc[wildcards.sample_name, "fastqdir"] + experimental_design.loc[wildcards.sample_name, "fastq1"],
        "read2": experimental_design.loc[wildcards.sample_name, "fastqdir"] + experimental_design.loc[wildcards.sample_name, "fastq2"]
    }
