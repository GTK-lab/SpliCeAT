#!/usr/bin/env python3

import sys
import gffutils
from snakemake.script import snakemake # type: ignore

sys.stderr = open(snakemake.log[0], "w")

gff = snakemake.input[0]
#gff = "resources/GRCm39_113_chr.gff3"
db = snakemake.output[0]
#db = "resources/GRCm39_119_chr.db"
gffutils.create_db(gff,db,
                force_gff=True,
                id_spec="ID",
                merge_strategy= "create_unique")
                #disable_infer_genes=True,
                #disable_infer_transcripts=True)