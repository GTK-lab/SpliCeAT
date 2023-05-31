rule get_novel_sequence:
    input:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf"
    output:
        novel=config["BASE_PATH"]+"results/augmented_transcriptome/merged_stringtie_assembly_novel_exon_filtered.fa",
        merged=config["BASE_PATH"]+"results/augmented_transcriptome/augmented_transcripts.fa"
    params:
        gffread=config["gff_read_path"],
        genome=config["genome_fasta"],
        transcripts=config["transcripts_fasta"]
    threads:
        4
    shell:
        "{params.gffread}" -w {output.novel} -g {params.genome} {input} && "
        "cat {params.transcripts} {output.novel} > {output.merged}"
