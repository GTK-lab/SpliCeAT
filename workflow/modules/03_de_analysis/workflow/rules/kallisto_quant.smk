rule kallisto_quant:
    input:
        unpack(get_fastq_files)
    output:
        "results/kallisto_quant_out/{sample_name}/abundance.tsv"
    params:
        INDEX = updated_file(aug_index,aug_output_dir),
        OUT_FILE = "results/kallisto_quant_out/{sample_name}",
        STRAND = lambda wildcards: "--rf-stranded" if config["experiment"]["strandedness"] == "reverse" else "--fr-stranded",
        BOOTSTRAPS = config["kallisto"]["bootstraps"],
    threads:
        config["kallisto"]["threads"]
    conda:
        "../../../../envs/kallisto.yaml"
    log:
        "logs/kallisto_quant/{sample_name}.log"
    message:
        "--- kallisto quant ---"
    shell:
        "kallisto quant -i {params.INDEX} -o {params.OUT_FILE} {params.STRAND} --bootstrap-samples={params.BOOTSTRAPS} --threads={params.N_THREADS} {input} > {log} 2>&1"
