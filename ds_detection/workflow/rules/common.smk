from snakemake.utils import validate
import pandas as pd

import yaml
from pathlib import Path

#validate(config, schema="../schemas/config.schema.yaml")

sample_file = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#")
annot = sample_file.set_index("sample_name")
annot['bam_dirs']= [str(f.parent) for f in [ Path(f) for f in annot['bam_file'] ]]
annot['bam_stem'] = [str(f.stem) for f in [ Path(f) for f in annot['bam_file'] ]]  # .stem removes the extension

leafcutter_grouppath = Path('results/leafcutter/groups.tsv')  
leafcutter_grouppath.parent.mkdir(parents=True, exist_ok=True)

sample_file.to_csv(leafcutter_grouppath,sep="\t",header=False,columns=["sample_name","group"],index=False)

# convert annot into dictionary for parametrization of rules, by deduplicating by sample_name (should only differ by bam_file)
samples = annot.to_dict(orient="index")

bam_files = list(annot['bam_file'])

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

def get_bam(wildcards):
    return annot.loc[wildcards.sample, "bam_file"]

def get_fq1(wildcards):
    return annot.loc[wildcards.sample, "fq1"]

def get_fq2(wildcards):
    return annot.loc[wildcards.sample, "fq2"]

def get_bam_stem(wildcards):
    return annot.loc[wildcards.sample, "bam_stem"]

def genome_release_name():
    build = config["ref"]["build"]
    release = config["ref"]["release"]
    flavor = config["ref"]["flavor"]
    if flavor:
        release += "_"
    return f"{build}_{release}{flavor}"

def gtf_file_path(filtered=False,gz=True):
    gz_str = ".gz" if gz else ""
    filt_str = "_filtered" if filtered else ""
    top_dir = "results" if filtered else "resources"    
    return f"{top_dir}/{genome_release_name()}{filt_str}.gtf{gz_str}"
    
def gff3_file_path(filtered=False,gz=True):
    gz_str = ".gz" if gz else ""
    filt_str = "_filtered" if filtered else ""
    top_dir = "results" if filtered else "resources"
    return f"{top_dir}/{genome_release_name()}{filt_str}.gff3{gz_str}"
    
def annotation_db_path():
    return f"resources/{genome_release_name()}.sqlite3"

def genome_file_path(gz=True):
    gz_str = ".gz" if gz else ""
    return f"resources/{genome_release_name()}.fa{gz_str}"
    
def leafcutter_output():
    if config['leafcutter']['activate']:
        return rules.leafcutter_differential_splicing.output

def majiq_output():
    if config['majiq']['activate']:
        return rules.majiq.output

import configparser

def majiq_files(group):
    stems = bamstems_for_group(group)
    files = [f"results/majiq/{f}.majiq" for f in stems]
    return files
    

def whippet_output():
    if config['whippet']['activate']:
        return rules.whippet.output

def leafcutter_strandedness(strandedness=config["experiment"]["strandedness"]):
    lookup = {
        "rf": "RF",
        "reverse": "RF",
        "forward": "FR",
        "fr": "FR",
        "unstranded": "XS",
        "yes": "FR"
    }
    return lookup[strandedness]
