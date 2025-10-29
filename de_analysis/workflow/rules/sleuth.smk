rule sleuth:
    input:
        expand(config["BASE_PATH"]+"/de_analysis/results/kallisto_quant_out/{sample_name}/abundance.tsv", sample_name = SAMPLES),
        config["BASE_PATH"]+"/augment_transcriptome/results/augmented_transcriptome/t2g_augment_uncollapsed.csv",
        config["BASE_PATH"]+"/augment_transcriptome/results/augmented_transcriptome/t2g_augment_collapsed.csv"
    output:
        config["BASE_PATH"]+"/de_analysis/results/sleuth/collapsed_differential_transcript_analysis_tpm.csv",
        config["BASE_PATH"]+"/de_analysis/results/sleuth/collapsed_differential_transcript_analysis.csv",
        config["BASE_PATH"]+"/de_analysis/results/sleuth/uncollapsed_differential_transcript_analysis_tpm.csv",
        config["BASE_PATH"]+"/de_analysis/results/sleuth/uncollapsed_differential_transcript_analysis.csv"
    params:
        script=config["BASE_PATH"]+"/de_analysis/scripts/sleuth.R"
    log:
        "logs/sleuth.log"
    conda:
        "../envs/sleuth.yaml"
    threads:
        12
    shell:
        "Rscript {params.script} > {log} 2>&1 "
