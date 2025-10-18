rule stringtie_merge:
    input:
        expand(config["BASE_PATH"]+"results/stringtie_assemblies/{sample}_ref_guided_assembly.gtf", sample=SAMPLES)
    output:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly.gtf"
    params:
        gtf=config["GTF"]
    threads:
        config["stringtie"]["threads"]
    conda:
        "../envs/stringtie.yaml"
    log:
        "logs/stringtie/stringtie_merge.log"
    shell:
        "stringtie --merge {input} -G {params.gtf} -i -o {output} > {log} 2>&1"
