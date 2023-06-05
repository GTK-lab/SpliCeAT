rule whippet_index:
    input:
        config["bam_dir"]+"{experiment}_treatment_merged.bam"
    output:
        config["BASE_PATH"]+"/results/index/{experiment}/{experiment}.jls"
    params:
        fasta=config["fasta_file"],
        gtf=config["annotation_gtf"],
        julia = config["julia"],
        whippet_bin = config["whippet_bin"]
    shell:
        "{params.julia} {params.whippet_bin}whippet-index.jl --fasta {params.fasta} --bam {input} --gtf {params.gtf} --bam-min-reads 20 --index {output}"
