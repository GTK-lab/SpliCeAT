import pandas as pd
from snakemake.script import snakemake # type: ignore
import logging
import sys

# Setup logging to Snakemake log file
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s'
)
# Load deltapsi TSV file
df = pd.read_csv(snakemake.input.dPSI_tsv, sep="\t")

# get strand information from het file
het_df = pd.read_csv(snakemake.input.het_tsv, sep="\t",comment='#')
strand_map = het_df[['lsv_id', 'strand','seqid']].drop_duplicates()

logging.info(f"dPSI LSV ID sample: '{df['lsv_id'].iloc[0]}'")
logging.info(f"HET LSV ID sample:  '{het_df['lsv_id'].iloc[0]}'")

# Merge strand information into the main dataframe
df = df.merge(strand_map, on='lsv_id', how='left')
df['strand'] = df['strand'].replace({1: '+', -1: '-', '1': '+', '-1': '-'})
df[['lsv_category', 'lsv_type_rest']] = df['lsv_type'].str.split('|', n=1, expand=True)

# Split 'lsv_type_rest' the same way as semicolon-separated columns
df['lsv_type_rest'] = df['lsv_type_rest'].astype(str).str.split('|')

# Identify columns that contain semicolon-separated values
vector_columns = [col for col in df.columns if df[col].astype(str).str.contains(';').any()]

# Split all vector columns into lists
for col in vector_columns:
	df[col] = df[col].astype(str).str.split(';')

vector_columns.append('lsv_type_rest')

# Ensure all list columns have the same length per row
# Explode the dataframe based on the vector columns
expanded_df = df.explode(vector_columns).reset_index(drop=True)
expanded_df['gene_id'] = expanded_df['gene_id'].str.removeprefix("gene:")
expanded_df['lsv_id'] = expanded_df['lsv_id'].str.removeprefix("gene:")

# Save the output to a new TSV file

expanded_df.to_csv(snakemake.output[0], sep="\t", index=False)

