#!/usr/bin/env python

import pandas as pd
from snakemake.script import snakemake # type: ignore


## First the effect sizes
df_effect = pd.read_csv(snakemake.input['effect_size'], sep="\t")

# Split the first column into four parts and extract the last character separately
df_effect[['chr', 'start', 'end', 'cluster']] = df_effect['intron'].str.rsplit(':', n=3, expand=True)
df_effect['strand'] = df_effect['cluster'].str[-1]  # Last character is the strand
df_effect['cluster'] = df_effect['cluster'].str[:-2]  # Remove last two characters (_+ or _-)

# Reorder columns for clarity
cols = ['chr', 'start', 'end', 'cluster', 'strand'] + \
    [col for col in df_effect.columns if col not in ['chr', 'start', 'end', 'cluster', 'strand', 'intron']]

df_effect = df_effect[cols]

## Now th clusters
# Load the TSV file
df_clusters = pd.read_csv(snakemake.input['signif'], sep="\t")

# Filter rows where status is "Success"
df_clusters = df_clusters[df_clusters['status'] == "Success"].copy()

# Split the 'cluster' column into three parts
cluster_split = df_clusters['cluster'].str.extract(r'(?P<chr>chr[0-9XYM]+):(?P<cluster>clu_\d+)_(?P<strand>[+-])')

# Concatenate the new columns with the existing dataframe
df_clusters = pd.concat([df_clusters.drop('cluster',axis=1), cluster_split], axis=1)

# Reorder columns for clarity
cols = ['chr', 'cluster', 'strand','loglr', 'df', 'p', 'p.adjust', 'genes']
df_clusters = df_clusters[cols]


lsvs = pd.merge(df_clusters,df_effect)
## reversing the insistence of leafcutter_ds.R
lsvs['chr'] = lsvs['chr'].str.removeprefix("chr")
cols = ['chr', 'start', 'end', 'cluster', 'strand'] + \
    [col for col in lsvs.columns if col not in ['chr', 'start', 'end', 'cluster', 'strand']]

lsvs = lsvs[cols]

lsvs.to_csv(snakemake.output['lsvs'],sep="\t",index=False)