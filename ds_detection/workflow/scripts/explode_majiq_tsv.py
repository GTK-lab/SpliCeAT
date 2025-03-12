import pandas as pd
from snakemake.script import snakemake # type: ignore

# Load the TSV file

df = pd.read_csv(snakemake.input[0], sep="\t")

# Identify columns that contain semicolon-separated values
vector_columns = [col for col in df.columns if df[col].astype(str).str.contains(';').any()]

# Split all vector columns into lists
for col in vector_columns:
    df[col] = df[col].astype(str).str.split(';')

# Ensure all list columns have the same length per row
# Explode the dataframe based on the vector columns
expanded_df = df.explode(vector_columns).reset_index(drop=True)

# Save the output to a new TSV file

expanded_df.to_csv(snakemake.output[0], sep="\t", index=False)

