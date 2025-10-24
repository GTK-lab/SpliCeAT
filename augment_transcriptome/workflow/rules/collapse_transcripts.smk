rule collapse_transcripts:
    input:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly.gtf",
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf"
    output:
        config["BASE_PATH"]+"results/augmented_transcriptome/t2g_augment_uncollapsed.csv",
        config["BASE_PATH"]+"results/augmented_transcriptome/t2g_augment_collapsed.csv"
    log:
        "logs/collapse_transcripts.log"
    conda:
       "../envs/for_R.yaml"
    params:
        script=config["BASE_PATH"]+"scripts/collapse_transcripts.R"
    threads:
        4
    shell:
        "Rscript {params.script} > {log} 2>&1"
