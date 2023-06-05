rule leafcutter_differential_splicing:
    input:
       counts=config["BASE_PATH"]+"results/{experiment}_perind_numers.counts.gz",
       pooled=config["BASE_PATH"]+"results/{experiment}_pooled",
       exons=config["BASE_PATH"]+"input/gtf_exons.txt.gz"
    output:
        config["BASE_PATH"]+"results/{experiment}_cluster_significance.txt",
        config["BASE_PATH"]+"results/{experiment}_effect_sizes.txt"
    params:
        groups_file = config["BASE_PATH"]+"config/{experiment}_groups_file.txt",
        prefix = config["BASE_PATH"]+"results/{experiment}",
        leafcutter_dir = config["leafcutter_dir"]
    threads:
        4
    shell:
        "{params.leafcutter_dir}scripts/leafcutter_ds.R --num_threads={threads} -e {input.exons} -i 1 -g 3 {input.counts} {params.groups_file} -o {params.prefix}"