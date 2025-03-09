

rule leafcutter:
    input:
        #leafcutter_ds="scripts/leafcutter_ds.R",
        #cluster_regtools="scripts/leafcutter_cluster_regtools.py",
        #exons = "results/leafcutter/exons.txt.gz",
        #cluster =  ["results/leafcutter/cluster_perind_numers.counts.gz",
        #            "results/leafcutter/cluster_pooled"],
        ds= ["results/leafcutter/cluster_significance.txt",
             "results/leafcutter/effect_sizes.txt"],
        pipeclean =  "results/leafcutter/.clean",
    output:
        "results/leafcutter/.donestamp",
    log:
        "logs/leafcutter/leafcutter.log"
    conda:
        "../envs/base.yaml"
    shell:
        "touch {output}"

rule cluster_regtools_script:
    input:
        storage.http("https://raw.githubusercontent.com/davidaknowles/leafcutter/refs/heads/master/clustering/leafcutter_cluster_regtools.py")
    output:
        "scripts/leafcutter_cluster_regtools.py"
    log:
        "logs/leafcutter/cluster_script.log"
    conda:
        "../envs/base.yaml"
    shell:
        "cp {input} {output}"


rule leafcutter_ds_script:
    input:
        storage.http("https://raw.githubusercontent.com/davidaknowles/leafcutter/refs/heads/master/scripts/leafcutter_ds.R")
    output:
        "scripts/leafcutter_ds.R"
    log:
        "logs/leafcutter_ds_script.log"
    conda:
        "../envs/base.yaml"
    shell:
        "cp {input} {output}"


rule leafcutter_exons:
    input:
        annotation_db_path(),
    params:
        biotype=config['experiment']['biotype'],
        max_support_level = config['experiment']['max_transcript_support_level'],
    output:
        "results/leafcutter/exons.txt.gz"
    conda:
        "../envs/gffutils.yaml"
    log:
        "logs/leafcutter/exons.log",
    script:
        "../scripts/leafcutter_exons.py"

rule regtools_junction_extract:
    input:
        lambda wildcards: get_bam(wildcards)
    output:
        "results/regtools/{sample}.junc"
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
        "results/leafcutter/junction_files.txt",
    log:
        "logs/leafcutter/junction_files.txt",
    run:
        wfile = open(str(output),"w")
        for line in input:
            wfile.write(f"{line}\n")


rule leafcutter_cluster:
    input:
        script = "scripts/leafcutter_cluster_regtools.py",
        junc_file =  "results/leafcutter/junction_files.txt",
    output:
        "results/leafcutter/cluster_perind_numers.counts.gz",
        "results/leafcutter/cluster_pooled",
        expand("{sample}.junc.cluster.sorted.gz",sample=samples),
    params:
        output = lambda w,output: os.path.basename(output[0]),
    log:
        "results/leafcutter/cluster.log"
    conda:
        "../envs/leafcutter.yaml",
    shell:
        "python3 {input.script} -j {input.junc_file} -o {params.output} -l 500000"
 

rule leafcutter_pipeclean:
    input:
        sorted=expand("{sample}.junc.cluster.sorted.gz",sample=samples),
        ds="results/leafcutter/cluster_significance.txt"
    output:
        "results/leafcutter/.clean"
    log:
        "logs/leafcutter/pipeclean.log",
    conda:
        "../envs/base.yaml"
    shell:
        "rm {input.sorted} && touch {output}"

rule leafcutter_differential_splicing:
    input:
        script="scripts/leafcutter_ds.R",
        counts="results/leafcutter/cluster_perind_numers.counts.gz",
        pooled="results/leafcutter/cluster_pooled",
        exons="results/leafcutter/exons.txt.gz",
    output:
        signif="results/leafcutter/cluster_significance.txt",
        effect_size="results/leafcutter/effect_sizes.txt"
    params:
        groups_file = leafcutter_grouppath,
        prefix = lambda w,output: os.path.dirname(output.effect_size),
    conda:
        "../envs/leafcutter.yaml"
    log:
        "logs/leafcutter/ds.log"
    threads:
        4
    shell:
        "Rscript {input.script} --num_threads={threads} -e {input.exons} -i 1 -g 3 {input.counts} {params.groups_file} -o {params.prefix} 2>{log} 1>&2;"
        "mv {params.prefix}_cluster_significance.txt {output.signif};"
        "mv {params.prefix}_effect_sizes.txt {output.effect_size}" 
