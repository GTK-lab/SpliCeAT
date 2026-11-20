#!/usr/bin/env python

import sys
import pandas as pd
# from snakemake.script import snakemake # type: ignore

sys.stderr = open(snakemake.log[0], "w")

g1, g2 = snakemake.params["groups"]
irf1, irf2 = snakemake.params["irf_groups"]
samples = snakemake.params["samples"]

diff = pd.read_csv(snakemake.input["all_results"], sep="\t", dtype={"ID": str})

# split the ID "GeneID;intron_type;chr:start-end" into columns (start is 0-based, BED style)
events = diff["ID"].str.extract(
	r"^(?P<gene_id>[^;]*);(?P<intron_type>[^;]*);(?P<chr>[^:]+):(?P<bed_start>\d+)-(?P<end>\d+)$"
)
events[["bed_start", "end"]] = events[["bed_start", "end"]].astype(int)

# Diff drops strand/gene name from its ids, this recover them from one sample's result
ir = pd.read_csv(snakemake.input["ir_result"], sep="\t", dtype={"Chr": str})
ir = pd.DataFrame({
	"chr": ir["Chr"],
	"bed_start": ir["Start"],
	"end": ir["End"],
	"gene_id": ir["Name"].str.split("/").str[1],
	"gene_name": ir["Name"].str.split("/").str[0],
	"strand": ir["Strand"],
}).drop_duplicates(subset=["chr", "bed_start", "end", "gene_id"])

events = events.merge(ir, on=["chr", "bed_start", "end", "gene_id"], how="left") 
print(f"{events['strand'].isna().sum()} of {len(events)} introns without a strand match", file=sys.stderr)
events["strand"] = events["strand"].fillna(".")

# ensure that log2fc is always grp 2 vs grp 1
if f"DESeq2.log2FoldChange.{irf2}_{irf1}" in diff.columns:		# already GROUPS[1] vs GROUPS[0]
	contrast = f"{irf2}_{irf1}"
	log2fc = diff[f"DESeq2.log2FoldChange.{contrast}"]
elif f"DESeq2.log2FoldChange.{irf1}_{irf2}" in diff.columns:    # GROUPS[0] vs GROUPS[1] -> flip
	contrast = f"{irf1}_{irf2}"
	log2fc = -diff[f"DESeq2.log2FoldChange.{contrast}"]
else:
	raise ValueError(f"[IRFinder]: no DESeq2 log2FoldChange column for {irf1}/{irf2}")

# build output table
out = pd.DataFrame({
	"chr": events["chr"],
	"start": events["bed_start"] + 1,  # BED (0-based) -> 1-based (ensure consistency with other tools)
	"end": events["end"],
	"strand": events["strand"],
	"gene_id": events["gene_id"].str.removeprefix("gene:").str.replace(r"\.[0-9]+$", "", regex=True),
	"gene_name": events["gene_name"],
	"intron_type": events["intron_type"],
	f"irratio_{g1}": diff[f"{irf1}.Mean.IRratio"],
	f"irratio_{g2}": diff[f"{irf2}.Mean.IRratio"],
	# change in IR ratio, group 2 - group 1 
	"delta_irratio": diff[f"{irf2}.Mean.IRratio"] - diff[f"{irf1}.Mean.IRratio"], # takes grp 2 - grp 1 
	"log2fc": log2fc,
	"padj": diff[f"DESeq2.padj.{contrast}"],
	"base_mean": diff[f"DESeq2.baseMean.{contrast}"],
})

# per-sample IR ratios
for i, sample in enumerate(samples, 1):
	out[f"irratio_{sample}"] = diff[f"IRratio.sample_{i}"]

out = out.sort_values(["chr", "start", "end"])
out.to_csv(snakemake.output["lsvs"], sep="\t", index=False)
print(f"wrote {len(out)} introns to {snakemake.output['lsvs']}", file=sys.stderr)


