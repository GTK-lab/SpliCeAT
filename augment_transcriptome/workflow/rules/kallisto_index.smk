rule kallisto_index:
    input:
        config["BASE_PATH"]+"results/augmented_transcriptome/augmented_transcripts.fa"
    output:
        config["BASE_PATH"]+"results/augmented_transcriptome/kallisto_index_augmented_transcriptome"
    threads:
        4
    shell:
        "kallisto index -i {output} {input}"
