rule leafcutter_cluster:
    input:
        unpack(get_output_junc)
    output:
        config["BASE_PATH"]+"/results/{experiment}_perind_numers.counts.gz",
        config["BASE_PATH"]+"/results/{experiment}_pooled"
    params:
        output = "results/{experiment}",
        out_dir = config["BASE_PATH"],
        junc_file = lambda wildcards: get_junc(wildcards),
        leafcutter_dir = config["leafcutter_dir"]
    shell:
        "python3 {params.leafcutter_dir}clustering/leafcutter_cluster_regtools.py -j {params.junc_file} -r {params.out_dir} -o {params.output} -l 500000"
