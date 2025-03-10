whippet_final_files = [
    "results/whippet/merged.bam",
    "results/whippet/merged.bam.bai",
    "results/whippet/whippet_index.jls",
    expand("results/whippet/quant/{sample}.psi.gz",sample=samples),
    "results/whippet/delta_psi.diff.gz",
]

comparison_groups = config["experiment"]["groups"]

rule whippet:
    input:
        whippet_final_files
    conda:
        "../envs/base.yaml"
    log:
        "logs/whippet.log"
    output:
        "results/whippet/.donestamp"
    shell:
        "touch {output}"


rule whippet_merge_bams:
    params:
        bam_groups = config['whippet']['bam_groups'],
    input:
        bam_files = bamfiles_for_groups(config['whippet']['bam_groups']),
    output:
        merged="results/whippet/merged.bam",
        merged_index="results/whippet/merged.bam.bai"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/whippet/merge_bams.log"
    threads: 8,
    shell:
        """
        samtools merge -o {output.merged} --threads {threads} {input.bam_files}; 
        samtools index {output.merged}
        """

rule whippet_index:
    input:
        genome_fa = genome_file_path(),
        bam =   "results/whippet/merged.bam",
        index = "results/whippet/merged.bam.bai",
        gtf = gtf_file_path(filtered=True),
    output:
        "results/whippet/whippet_index.jls"
    params:
        bam_min_reads = config["whippet"]["bam_min_reads"]
    log:
        "logs/whippet/whippet_index.log"
    conda:
        "../envs/whippet.yaml"
    shell:
        "julia ${{WHIPPET_PATH}}/whippet-index.jl "
        "--fasta {input.genome_fa} "
        "--bam {input.bam} "
        "--gtf {input.gtf} "
        "--bam-min-reads {params.bam_min_reads} "
        "--index {output} >{log} 2>&1;"
        
rule whippet_quantify:
    input:
        fq1 = get_fq1,
        fq2 = get_fq2,
        index = "results/whippet/whippet_index.jls",
    output:
        "results/whippet/quant/{sample}.psi.gz",
    params:
        out = lambda w, output: output[0].removesuffix(".psi.gz")
    log:
        "logs/whippet/{sample}_quant.log",
    conda:
        "../envs/whippet.yaml",
    shell:
        "julia ${{WHIPPET_PATH}}/whippet-quant.jl -o {params.out} -x {input.index} {input.fq1} {input.fq2} 2>{log}"


rule whippet_delta:
    input:
        grp1=expand("results/whippet/quant/{sample}.psi.gz",sample=sample_for_group(config["experiment"]["groups"][0])),
        grp2=expand("results/whippet/quant/{sample}.psi.gz",sample=sample_for_group(config["experiment"]["groups"][1]))
    params:
        grp1=",".join(expand("results/whippet/quant/{sample}.psi.gz",sample=sample_for_group(config["experiment"]["groups"][0]))),
        grp2=",".join(expand("results/whippet/quant/{sample}.psi.gz",sample=sample_for_group(config["experiment"]["groups"][1]))),
    output:
        "results/whippet/delta_psi.diff.gz",
    log:
        "logs/whippet/delta_psi.log",
    conda:
        "../envs/whippet.yaml",
    shell:
        "julia ${{WHIPPET_PATH}}/whippet-delta.jl "
        "-a {params.grp1} -b {params.grp2} -o results/whippet/delta_psi "
