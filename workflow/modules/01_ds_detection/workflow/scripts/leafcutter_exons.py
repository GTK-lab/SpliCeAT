#!/usr/bin/env python3
import sys
import pybedtools
import gzip
from snakemake.script import snakemake # type: ignore

sys.stderr = open(snakemake.log[0], "w")

def strip_prefix(feature):
    feature.name = feature.name.removeprefix("gene:")
    return feature

gtf = pybedtools.BedTool(snakemake.input[0])
exons = gtf.filter(lambda x: x[2] == 'exon')

if snakemake.output[0].endswith(".gz"):
    sys.stdout = gzip.open(snakemake.output[0],"wt")
else:
    sys.stdout = open(snakemake.output[0],"w")

print(*['chr','start', 'end', 'strand', 'gene_name'],sep='\t')

for exon in exons:
    exon.name = exon.attrs['gene_id'].removeprefix("gene:")
    print(*[exon.chrom,exon.start+1,exon.stop,exon.strand,exon.name],sep='\t')