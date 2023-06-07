rule stringtie_merge:
    input:
        unpack(get_stringtie_assembly)
    output:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly.gtf"
    params:
        stringtie=config["STRINGTIE_COMMAND"],
        gtf=config["GTF"]
    threads:
        config["STRINGTIE_THREADS"]
    shell:
        "{params.stringtie} --merge {input} -G {params.gtf} -i -o {output}"
 
rule stringtie_merge:
    input:
        config["BASE_PATH"]+"results/stringtie_assemblies/{sample}_ref_guided_assembly.gtf"
    output:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly.gtf"
    params:
        stringtie=config["STRINGTIE_COMMAND"],
        gtf=config["GTF"]
    threads:
        config["STRINGTIE_THREADS"]
    shell:
        "{params.stringtie} --merge {input} -G {params.gtf} -i -o {output}"
 
