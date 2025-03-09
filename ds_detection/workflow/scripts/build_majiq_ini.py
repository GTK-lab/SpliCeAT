#!/usr/bin/env python

import configparser
from snakemake.script import snakemake # type: ignore
import sys
import pandas as pd
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")
output_file = open(snakemake.output[0],"w")

samples = pd.read_csv(snakemake.config["samples"],sep="\t", dtype=str, comment="#")
bam_files = [ Path(f) for f in samples['bam_file'] ]
bam_dirs = list(set([str(f.parent) for f in bam_files]))
bam_dirs = [ f + "/" for f in bam_dirs ]
print(bam_dirs)

samples['bam_name'] = [str(f.stem) for f in bam_files]  # .stem removes the extension
grouped_bam_names = samples.groupby("group")["bam_name"].agg(",".join).to_dict()

majiq_config = configparser.ConfigParser()

majiq_config['info'] = {
        "experiment" : snakemake.config["experiment"]["name"],
        "genome" : snakemake.config["ref"]["build"],
        "bamdirs" : ",".join(bam_dirs),
        "strandedness" : snakemake.config["experiment"]["strandedness"]
    }

majiq_config['experiments'] = grouped_bam_names

majiq_config.write(output_file)