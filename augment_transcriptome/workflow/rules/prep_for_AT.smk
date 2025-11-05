rule prep_for_AT:
    input:
        config["majiq_lsv_file_path"],
        config["leafcutter_cluster_sig_file_path"],
        config["leafcutter_effect_size_file_path"],
        config["whippet_diff_file_path"],
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly.gtf"
    output:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf"
    log:
        "logs/prep_for_AT.log"
    conda:
        "../envs/for_R.yaml"
    params:
        script=config["BASE_PATH"]+"scripts/prep_for_AT.R"
    threads:
        4
    shell:
        "Rscript {params.script} > {log} 2>&1"
