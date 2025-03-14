
rule get_genome:
    output:
        temp(genome_file_path(gz=False)),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        datatype="dna",
    cache: True,
    log:
        "logs/get_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v5.8.2/bio/reference/ensembl-sequence"

rule gzip_genome:
    input:
        genome_file_path(gz=False),
    output:
        genome_file_path()
    threads: 1
    log:
        "logs/bgzip/genome.log",
    wrapper:
        "v5.8.2/bio/bgzip"        


rule get_annotation_gff:
    output:
        gff3_file_path(),
        #"resources/annotation.gtf.gz",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="chr"
    log:
        "logs/get_annotation.log",
    cache: True #"omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v5.8.2/bio/reference/ensembl-annotation"

rule get_annotation_gtf:
    output:
        gtf_file_path(),
        #"resources/annotation.gtf.gz",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="chr",
    log:
        "logs/get_annotation-gtf.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v5.8.2/bio/reference/ensembl-annotation"


rule make_annotation_db:
    input:
        gff3_file_path(),
    output:
        annotation_db_path(),
    log:
        "logs/annotation_db.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/make_annotation_db.py"


rule make_filtered_gff:
    input:
        annotation_db_path(),
    params:
        biotype=config['experiment']['biotype'],
        max_support_level = config['experiment']['max_transcript_support_level'],
    output:
        gff3_file_path(filtered=True,gz=False)
    conda:
        "../envs/gffutils.yaml"
    log:
        "logs/ref/filter_gff.log"
    script:
        "../scripts/filter_annotation_to_gff.py"
    

rule make_filtered_gtf:
    input:
        gff3_file_path(filtered=True,gz=False),
    output:
        gtf_file_path(filtered=True),
    conda:
        "../envs/gffread.yaml",
    log:
        "logs/ref/filter_gtf.log"
    shell:
        "gffread {input} -T | bgzip > {output} 2>{log}"

rule uncompress:
    wildcard_constraints:
        filename=".*(?<!\\.gz)"
    input:
        "{filename}.gz"
    log:
        "logs/uncompress_{filename}.log"
    output:
        "{filename}",
    wrapper:
        "v5.8.3/bio/bgzip"