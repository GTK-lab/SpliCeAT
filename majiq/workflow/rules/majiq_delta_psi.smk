rule majiq_delta_psi:
    input:
        config["BASE_PATH"]+"/results/majiq_build/{experiment}/majiq.log"
    output:
        config["BASE_PATH"]+"/results/majiq_delta_psi/{experiment}/ctr-treatment.deltapsi.tsv"
    params:
        out=config["delta_out"],
        grp1=lambda wildcards: delta_psi_grp1(wildcards),
        grp2=lambda wildcards: delta_psi_grp2(wildcards)
    threads:
        16
    shell:
        "majiq deltapsi -grp1 {params.grp1} -grp2 {params.grp2} -j {threads} -o {params.out} -n ctr treatment"
