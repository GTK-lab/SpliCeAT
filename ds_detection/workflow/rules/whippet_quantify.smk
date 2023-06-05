rule whippet_quantify:
    input:
        fq1 = lambda wildcards: get_fastq(wildcards)[0],
        fq2 = lambda wildcards: get_fastq(wildcards)[1],
        index = lambda wildcards: get_index(wildcards)
    output:
        config["BASE_PATH"]+"/results/quantify/{sample}/{sample}.psi.gz"
    params:
        out = config["BASE_PATH"]+"/results/quantify/{sample}/{sample}",
        julia = config["julia"],
        whippet_bin = config["whippet_bin"]
    shell:
        "{params.julia} {params.whippet_bin}whippet-quant.jl {input.fq1} {input.fq2} -o {params.out} -x {input.index}"
