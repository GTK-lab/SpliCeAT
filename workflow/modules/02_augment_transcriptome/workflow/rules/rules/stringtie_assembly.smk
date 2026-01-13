rule stringtie_assembly:
    input:
        unpack(get_bam_files)
    output:
        config["BASE_PATH"]+"results/stringtie_assemblies/{sample}_ref_guided_assembly.gtf"
    params:
        gtf=config["GTF"],
        juncs=config["stringtie"]["junc_cutoff"],
        threads=config["stringtie"]["threads"],
        strand=config["stringtie"]["strand"]
    conda:
        "../envs/stringtie.yaml"
    log:
        "logs/stringtie/{sample}_stringtie_assembly.log"
    threads:
        config["stringtie"]["threads"]
    shell:
        "stringtie {input} -G {params.gtf} -o {output} -j {params.juncs} -p {params.threads} --{params.strand} > {log} 2>&1"
