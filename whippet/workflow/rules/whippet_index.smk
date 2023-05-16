rule whippet_index:
    input:
        config["bam_dir"]+"{experiment}_treatment_merged.bam"
    output:
        config["base_path"]+"/results/index/{experiment}/{experiment}.jls"
    params:
        fasta=config["fasta_file"],
        gtf=config["annotation_gtf"]
    shell:
        "~/julia-1.7.2/bin/julia ~/Whippet.jl/bin/whippet-index.jl --fasta {params.fasta} --bam {input} --gtf {params.gtf} --bam-min-reads 20 --index {output}"