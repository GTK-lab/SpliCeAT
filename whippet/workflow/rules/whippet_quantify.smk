rule whippet_quantify:
    input:
        fq1 = lambda wildcards: get_fastq(wildcards)[0],
        fq2 = lambda wildcards: get_fastq(wildcards)[1],
        index = lambda wildcards: get_index(wildcards)
    output:
        config["base_path"]+"/results/quantify/{sample}/{sample}.psi.gz"
    params:
        config["base_path"]+"/results/quantify/{sample}/{sample}"
    shell:
        "~/julia-1.7.2/bin/julia ~/Whippet.jl/bin/whippet-quant.jl {input.fq1} {input.fq2} -o {params} -x {input.index}"