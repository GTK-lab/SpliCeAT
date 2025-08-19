

rule leafcutter_exons:
    input:
        gtf_file_path(filtered=True)
    output:
        "results/leafcutter/exons.txt.gz",
    conda:
        "../envs/gffutils.yaml"
    log:
        "logs/leafcutter/exons.log",
    script:
        "../scripts/leafcutter_exons.py"


rule regtools_junction_extract:
    input:
        get_bam,
    output:
        "results/regtools/{sample}.junc",
    container:
        "docker://griffithlab/regtools:release-1.0.0"
    conda:
        "../envs/regtools.yaml"
    log:
        "logs/regtools/{sample}.log"
    params:
        strand = leafcutter_strandedness(),
    shell:
        "echo strand: {params.strand} 2>{log} 1>&2;"
        "regtools junctions extract -a 8 -m 50 -M 500000 -s {params.strand} {input} -o {output} 2>>{log}"


rule create_leafcutter_junction_file:
    input:
        expand("results/regtools/{sample}.junc",sample=samples),
    output:
        temp("results/leafcutter/junction_files.txt"),
    log:
        "logs/leafcutter/junction_files.txt",
    run:
        wfile = open(str(output),"w")
        for line in input:
            wfile.write(f"{line}\n")


rule leafcutter_cluster:
    input:
        rules.create_leafcutter_junction_file.output,
    output:
        counts="results/leafcutter/leafcutter_perind_numers.counts.gz",
        pooled="results/leafcutter/leafcutter_pooled",
    #shadow: "full",
    params:
        max_intron_length = config["leafcutter"]["max_intron_length"],
        rundir=lambda w,output: os.path.dirname(output.counts),
    log:
        "logs/leafcutter/cluster.log"
    conda:
        "../envs/leafcutter.yaml",
    shell:
        "leafcutter_cluster_regtools.py -j {input} -r {params.rundir} -l {params.max_intron_length}"
 


rule leafcutter_differential_splicing:
    input:
        counts="results/leafcutter/leafcutter_perind_numers.counts.gz",
        pooled="results/leafcutter/leafcutter_pooled",
        exons="results/leafcutter/exons.txt.gz",
    output:
        signif="results/leafcutter/cluster_significance.txt",
        effect_size="results/leafcutter/effect_sizes.txt"
    params:,
        min_samples_per_group=config["leafcutter"]["min_samples_per_group"],
        groups_file = leafcutter_grouppath,
        prefix = lambda w,output: os.path.dirname(output.effect_size),
    conda:
        "../envs/leafcutter.yaml"
    log:
        "logs/leafcutter/ds.log"
    threads:
        4
    shell:
        "unset R_LIBS_USER; "
        "leafcutter_ds.R --num_threads={threads} -e {input.exons} -i 1 -g {params.min_samples_per_group} {input.counts} {params.groups_file} -o {params.prefix} 2>{log} 1>&2;"
        "mv {params.prefix}_cluster_significance.txt {output.signif};"
        "mv {params.prefix}_effect_sizes.txt {output.effect_size}" 


rule leafcutter_lsvs:
    input:
        signif="results/leafcutter/cluster_significance.txt",
        effect_size="results/leafcutter/effect_sizes.txt",
    output:
        lsvs="results/leafcutter/leafcutter_lsvs.tsv",
    log:
        "logs/leafcutter/lsvs.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/leafcutter_lsvs.py"