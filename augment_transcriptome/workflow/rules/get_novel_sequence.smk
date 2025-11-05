rule get_novel_sequence:
    input:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf"
    output:
        novel=config["BASE_PATH"]+"results/augmented_transcriptome/merged_stringtie_assembly_novel_exon_filtered.fa",
        merged=config["BASE_PATH"]+"results/augmented_transcriptome/augmented_transcripts.fa"
    params:
        genome=config["genome_fasta"],
        transcripts=config["transcripts_fasta"]
    conda:
        "../../../ds_detection/workflow/envs/gffread.yaml"
    threads:
        4
    shell:
        "gffread -w {output.novel} -g {params.genome} {input} && "
        "cat {params.transcripts} {output.novel} > {output.merged}"
