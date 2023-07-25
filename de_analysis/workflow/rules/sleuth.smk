rule sleuth:
    input:
        expand(config["BASE_PATH"]+"/de_analysis/results/kallisto_quant_out/{sample_name}/abundance.tsv", sample_name = SAMPLES)
    output:
        config["BASE_PATH"]+"/de_analysis/results/sleuth/collapsed_differential_transcript_analysis_tpm.csv",
        config["BASE_PATH"]+"/de_analysis/results/sleuth/collapsed_differential_transcript_analysis.csv",
        config["BASE_PATH"]+"/de_analysis/results/sleuth/uncollapsed_differential_transcript_analysis_tpm.csv",
        config["BASE_PATH"]+"/de_analysis/results/sleuth/uncollapsed_differential_transcript_analysis.csv"
    params:
        script=config["BASE_PATH"]+"/de_analysis/scripts/sleuth.R"
    threads:
        12
    shell:
        "Rscript {params.script}"
