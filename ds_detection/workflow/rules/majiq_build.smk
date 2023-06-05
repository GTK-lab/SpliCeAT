rule majiq_build:
    input:
        unpack(get_conf_file)
    output:
        config["BASE_PATH"]+"/results/majiq_build/{experiment}/majiq.log"
    params:
        gtf=config["gtf"],
        out=config["build_out"]
    threads:
        8
    shell:
        "majiq build {params.gtf} -c {input} -j {threads} -o {params.out}"
