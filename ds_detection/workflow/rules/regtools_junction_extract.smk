rule regtools_junction_extract:
    input:
        lambda wildcards: get_bam(wildcards)
    output:
        config["BASE_PATH"]+"results/{sample}.junc"
    params:
        strand = config["STRAND"]
    shell:
        "regtools junctions extract -s {params.strand} {input} -o {output}"