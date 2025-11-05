#!/usr/bin/env python3
import sys
import gffutils
import gzip
from snakemake.script import snakemake # type: ignore

sys.stderr = open(snakemake.log[0], "w")
transcript_type = snakemake.config["experiment"]["biotype"]

db = gffutils.FeatureDB(snakemake.input[0],keep_order=True)

if snakemake.output[0].endswith(".gz"):
    output_file = gzip.open(snakemake.output[0],"wt")
else:
    output_file = open(snakemake.output[0],"w")

with output_file as output:
    for f in db.features_of_type('exon'):
        if (not transcript_type) or (transcript_type  == f['transcript_biotype'][0]):
            print(*[f.seqid,f.start,f.end,f['gene_id'][0],1000,f.strand],sep='\t',file=output)

