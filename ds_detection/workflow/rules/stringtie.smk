rule stringtie_assembly:
    input:
        bam = get_bam,
        gff=gff3_file_path(filtered=True,gz=False),
    output:
        "results/stringtie/assemblies/{sample}.gtf"
    params:
        juncs=config["stringtie"]["min_juncs"],
        strand=config["stringtie"]["strand"]
    threads:
        config["stringtie"]["threads"]
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie {input.bam} -G {input.gff} -o {output} -j {params.juncs} -p {threads} --{params.strand}"    


rule stringtie_merge:
    input:
        assemblies=expand("results/stringtie/assemblies/{sample}.gtf",sample=samples),
        gff=gff3_file_path(filtered=True,gz=False)        
    output:
        "results/stringtie/merged.gtf"
    threads:
        config["stringtie"]["threads"],
    conda:
        "../envs/stringtie.yaml"
    shell:
        "stringtie --merge {input.assemblies} -G {input.gff} -i -o {output}"