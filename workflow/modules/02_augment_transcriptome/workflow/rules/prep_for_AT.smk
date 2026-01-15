#transcriptome augmentation
rule prep_for_AT:
    input:
        ds_files = copied_outputs(ds_output, ds_output_dir),
        gtf_merged="results/merged_assembly/merged_stringtie_assembly.gtf"
    output:
        merged_filtered_gtf="results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf",
        merged_fil_withRef_gtf="results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered_with_reference.gtf"
    log:
        "logs/prep_for_AT.log"
    conda:
        "../../../../envs/for_R.yaml"
    params:
        module_path = os.path.join(config["pipeline_dir"], "workflow/modules/02_augment_transcriptome/"),
        organism = config["ref"]["species"],
        ensembl=config["ref"]["ensembl_release"]
    threads:
        4
    script:
        "../scripts/prep_for_AT.R"

rule collapse_transcripts:
    input:
        gtf_merged="results/merged_assembly/merged_stringtie_assembly.gtf",
        gtf_filtered="results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf"
    output:
        uncollapsed="results/augmented_transcriptome/t2g_augment_uncollapsed.csv",
        collapsed="results/augmented_transcriptome/t2g_augment_collapsed.csv"
    log:
        "logs/collapse_transcripts.log"
    conda:
       "../../../../envs/for_R.yaml"
    params:
        module_path = os.path.join(config["pipeline_dir"], "workflow/modules/02_augment_transcriptome/"),
        organism = config["ref"]["species"],
        ensembl=config["ref"]["ensembl_release"]
    threads:
        4
    script:
        "../scripts/collapse_transcripts.R"

rule get_novel_sequence:
    input:
        "results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf"
    output:
        novel="results/augmented_transcriptome/merged_stringtie_assembly_novel_exon_filtered.fa",
        merged="results/augmented_transcriptome/augmented_transcripts.fa"
    params:
        genome=updated_file(genome_file_path(gz=False,dna=True),ref_target_dir),
        transcripts=updated_file(genome_file_path(gz=False,dna=False),ref_target_dir)
    conda:
        "../../../../envs/genomic_utils.yaml"
    threads:
        4
    shell:
        "gffread -w {output.novel} -g {params.genome} {input} && "
        "cat {params.transcripts} {output.novel} > {output.merged}"