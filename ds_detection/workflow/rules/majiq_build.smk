rule majiq_build:
    input:
        unpack(get_conf_file)
    output:
        config["BASE_PATH"]+"/results/majiq_build/{experiment}/majiq.log"
    params:
        gff3=config["gff3"],
        out=config["build_out"]
    threads:
        8
    shell:
        "majiq build {params.gff3} -c {input} -j {threads} -o {params.out}"
