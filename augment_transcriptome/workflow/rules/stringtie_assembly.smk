rule stringtie_assembly:
    input:
        unpack(get_bam_files)
    output:
        config["BASE_PATH"]+"results/stringtie_assemblies/{sample}_ref_guided_assembly.gtf"
    params:
        stringtie=config["STRINGTIE_COMMAND"],
        gtf=config["GTF"],
        juncs=config["JUNCS_CUTOFF"],
        threads=config["STRINGTIE_THREADS"],
        strand=config["STRAND"]
    shell:
        "{params.stringtie} {input} -G {params.gtf} -o {output} -j {params.juncs} -p {params.threads} --{params.strand}"    