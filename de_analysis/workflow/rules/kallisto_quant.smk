rule kallisto_quant:
    input:
        unpack(get_fastq_files)
    output:
        config["OUT_FILE"]+"/abundance.tsv"
    params:
        INDEX = config["INDEX"],
        OUT_FILE = config["OUT_FILE"],
        STRAND = config["STRAND"],
        BOOTSTRAPS = config["BOOTSTRAPS"],
        N_THREADS = config["N_THREADS"]
    threads:
        config["N_THREADS"]
    message:
        "--- kallisto quant ---"
    shell:
        "kallisto quant -i {params.INDEX} -o {params.OUT_FILE} --{params.STRAND} --bootstrap-samples={params.BOOTSTRAPS} --threads={params.N_THREADS} {input}"