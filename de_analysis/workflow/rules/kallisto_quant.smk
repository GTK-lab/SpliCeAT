rule kallisto_quant:
    input:
        unpack(get_fastq_files)
    output:
        config["BASE_PATH"]+"/de_analysis/results/kallisto_quant_out/{sample_name}/abundance.tsv"
    params:
        INDEX = config["BASE_PATH"]+"/augment_transcriptome/results/augmented_transcriptome/kallisto_index_augmented_transcriptome",
        OUT_FILE = config["BASE_PATH"]+"/de_analysis/results/kallisto_quant_out/{sample_name}",
        STRAND = config["STRAND"],
        BOOTSTRAPS = config["BOOTSTRAPS"],
        N_THREADS = config["N_THREADS"]
    threads:
        config["N_THREADS"]
    conda:
        "../../../augment_transcriptome/workflow/envs/kallisto.yaml"
    log:
        "logs/kallisto_quant/{sample_name}.log"
    message:
        "--- kallisto quant ---"
    shell:
        "kallisto quant -i {params.INDEX} -o {params.OUT_FILE} --{params.STRAND} --bootstrap-samples={params.BOOTSTRAPS} --threads={params.N_THREADS} {input} > {log} 2>&1"
