
majiq_sj_files = [ f"results/majiq/{f}.sj" for f in annot['bam_stem'] ]
majiq_majiq_files = [ f"results/majiq/{f}.majiq" for f in annot['bam_stem'] ]

comparison_groups = config["experiment"]["groups"]

majiq_final_files = [
    expand("results/majiq/{samples.bam_stem}.sj",samples=annot.itertuples),
    f"results/majiq/{'-'.join(comparison_groups)}.deltapsi.tsv",
    f"results/majiq/{'-'.join(comparison_groups)}.deltapsi.voila",
    expand("results/majiq/{samples.bam_stem}.majiq",samples=annot.itertuples)
]

rule majiq:
    input:
        majiq_final_files
    output:
        "results/majiq/.donestamp"
    log:
        "logs/majiq.log"
    conda:
        "../envs/base.yaml"
    shell:
        "touch {output}"

rule majiq_conf:
    log:
        "logs/majiq/build_ini.log",
    output:
        "config/majiq/majiq.ini"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/build_majiq_ini.py"


rule majiq_build:
    input:
        gff3=gff3_file_path(filtered=True),
        conf_file = "config/majiq/majiq.ini",
        bams = bam_files,
    params:
        license=config["majiq"]["license"],
    output:
        "results/majiq/splicegraph.sql",
        majiq_sj_files,
        majiq_majiq_files
    log:
        "logs/majiq/majiq_build.log",
    conda:
        "../envs/majiq.yaml"
    threads:
        8
    shell:
        "majiq --license {params.license} build {input.gff3} "
        "--logger {log} "
        "-c {input.conf_file} -j {threads} -o results/majiq >/dev/null 2>&1" 


rule majiq_delta_psi:
    input:
        sj_files=majiq_sj_files,
        majiq_files=majiq_majiq_files,
        splicegraph="results/majiq/splicegraph.sql",
    output:
        tsv=f"results/majiq/{'-'.join(comparison_groups)}.deltapsi.tsv",
        voila=f"results/majiq/{'-'.join(comparison_groups)}.deltapsi.voila",
    log:
        "logs/majiq/majiq_deltapsi.log",
    params:
        output_dir=lambda w,input: os.path.dirname(input.splicegraph),
        license=config['majiq']['license'] ,
        names=lambda _: " ".join(comparison_groups),
        prefix=lambda _: "-".join(comparison_groups),
        grp1= lambda _: " ".join(majiq_files(comparison_groups[0])),
        grp2= lambda _: " ".join(majiq_files(comparison_groups[1])),
    threads:
        16
    conda:
        "../envs/majiq.yaml"
    shell:
        "majiq --license {params.license} deltapsi "
        "-grp1 {params.grp1} -grp2 {params.grp2} "
        "--logger {log} "
        "-j {threads} -o {params.output_dir} -n {params.names}"
