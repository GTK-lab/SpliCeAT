from snakemake.utils import validate
import pandas as pd
import os
import yaml
import configparser
from pathlib import Path


# forcing leafcutter group order to match config file
group_order = config["experiment"]["groups"]
sample_file_df['group'] = pd.Categorical(
    sample_file_df['group'],
    categories=group_order,
    ordered=True
)
sample_file_df = sample_file_df.sort_values("group")

leafcutter_grouppath = Path(os.path.join(LC_DIR,'leafcutter_groups.tsv'))
leafcutter_grouppath.parent.mkdir(parents=True, exist_ok=True)

sample_file_df.to_csv(
    leafcutter_grouppath,
    sep="\t",
    header=False,
    columns=["sample_name", "group"],
    index=False
)

def leafcutter_output():
	if config['leafcutter']['activate']:
		return rules.leafcutter_lsvs.output

def majiq_output():
	if config['majiq']['activate']:
		return rules.majiq.output

def majiq_files(group):
	stems = bamstems_for_group(group)
	files = [os.path.join(R1_DS,f"majiq/{f}.majiq") for f in stems]
	return files

def whippet_output():
	if config['whippet']['activate']:
		return rules.whippet.input

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
