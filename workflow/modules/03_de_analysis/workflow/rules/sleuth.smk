rule sleuth:
    input:
        aug_files=copied_outputs(aug_output,aug_output_dir),
        kallisto_tsv=expand("results/kallisto_quant_out/{sample_name}/abundance.tsv", sample_name = SAMPLES)
    output:
        result_collapsed_tpm="results/sleuth/collapsed_differential_transcript_analysis_tpm.csv",
        result_collapsed="results/sleuth/collapsed_differential_transcript_analysis.csv",
        result_uncollapsed_tpm="results/sleuth/uncollapsed_differential_transcript_analysis_tpm.csv",
        result_uncollapsed="results/sleuth/uncollapsed_differential_transcript_analysis.csv"
    params:
        module_path= os.path.join(config["pipeline_dir"], "workflow/modules/03_de_analysis/"),
        design_matrix= os.path.join(config["pipeline_dir"], "config/",config["samples"])
    log:
        "logs/sleuth.log"
    conda:
        "../../../../envs/sleuth.yaml"
    threads:
        12
    script:
        "../scripts/sleuth.R"