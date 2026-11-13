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

def irfinder_output():
	if config['irfinder']['activate']:
		return rules.irfinder_lsvs.output 

def leafcutter_strandedness(strandedness=config["experiment"]["strandedness"]):
	val = str(strandedness).lower()
	lookup = {
		"yes": "FR",
		"forward": "FR",
		"fr": "FR",
		"rf": "RF",
		"reverse": "RF",
		"none" : "XS",
		"no" : "XS",
		"unstranded": "XS"
	}
	return lookup.get(val,"XS")

# check sample names and make group names safe for IRFinder Diff
def prepare_irfinder_names(groups, sample_for_group):

    # Check sample names
    for sample in sample_for_group(groups[0]) + sample_for_group(groups[1]):
        if "-h" in sample:
            raise ValueError(
                f"[IRFinder]: sample name '{sample}' contains '-h', "
                "which IRFinder Diff reads as its help option."
            )

    # Check group names 
    irf_groups = []
    for group in groups[:2]:
        name = str(group).lower()
        safe = "".join(
            c if (c.isalnum() or c == ".") else "."
            for c in name
        )
        if not safe or not safe[0].isalpha():
            safe = "g." + safe
        irf_groups.append(safe)

    # Check that groups are not identical (to avoid being treated as the same group by IRFinder Diff)
    if irf_groups[0] == irf_groups[1]:
        raise ValueError(
            f"[IRFinder]: groups '{groups[0]}' and '{groups[1]}' "
            f"both become '{irf_groups[0]}'."
        )

    return irf_groups

# build list of IRFinder result files for a given group
def irfinder_group_files(group):
	return " ".join(os.path.join("quant", s, IRF_RESULT) for s in sample_for_group(group))
