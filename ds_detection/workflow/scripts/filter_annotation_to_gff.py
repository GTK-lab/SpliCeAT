#!/usr/bin/env python3

import sys
import gffutils
import gzip
from snakemake.script import snakemake # type: ignore

sys.stderr = open(snakemake.log[0], "w")
biotype = snakemake.params["biotype"]
max_transcript_support = snakemake.params["max_support_level"]

db = gffutils.FeatureDB(snakemake.input[0],keep_order=True)

if snakemake.output[0].endswith(".gz"):
    fout = gzip.open(snakemake.output[0],"wt")
else:
    fout = open(snakemake.output[0],"w")


for d in db.directives:
    fout.write("#" + str(d) + '\n')

for g in db.features_of_type('gene'):
    if (not biotype) or (biotype  == g['biotype'][0]):
        gid = g.id.removeprefix("gene:")
        fout.write(str(g) + '\n')
        for tx in db.children(g,featuretype=('transcript','mRNA')):
            if (not biotype) or (biotype  == tx['biotype'][0]):
                if tx.attributes.get('transcript_support_level'):
                    TSL = tx.attributes['transcript_support_level'][0].split(' ', 1)[0]
                    if TSL.isdigit() and int(TSL) <= max_transcript_support:
                        fout.write(str(tx) + '\n')
                        for exon in db.children(tx,featuretype="exon"):
                            fout.write(str(exon) + '\n')

