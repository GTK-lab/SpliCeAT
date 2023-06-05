rule create_genome_txt:
    input:
       gtf=config["annotation_gtf"]
    output:
        config["BASE_PATH"]+"/input/gtf_exons.txt.gz"
    params:
        leafcutter_dir = config["leafcutter_dir"]
    threads:
        4
    shell:
        "{params.leafcutter_dir}scripts/gtf_to_exons.R {input} {output}"
