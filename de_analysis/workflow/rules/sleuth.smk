rule sleuth:
    input:
        config["BASE_PATH"]+"/de_analysis/results/kallisto_quant_out/{sample_name}/abundance.tsv"
    output:
        config["BASE_PATH"]+"/de_analysis/results/kallisto_quant_out/differential_transcript_expression_analysis_uncollapsed.csv",
        config["BASE_PATH"]+"/de_analysis/results/kallisto_quant_out/differential_transcript_expression_analysis_collapsed.csv"
    params:
        script=config["BASE_PATH"]+"/de_analysis/scripts/sleuth.R"
    threads:
        4
    shell:
        "Rscript {params.script}"
