rule get_novel_sequence:
    input:
        config["BASE_PATH"]+"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf"
    output:
        config["BASE_PATH"]+"results/augmented_transcriptome/merged_stringtie_assembly_novel_exon_filtered.fa"
    params:
        gffread=config["gff_read_path"],
        genome=config["genome_fasta"]
    threads:
        4
    shell:
        "{params.gffread}"
        
        
    ~/gffread/gffread \
-w /mnt/cbis/home/yongshan/at_pipeline/temp_results/stringtie_filtering/exon_level/nestin_ctx_e14_ref_guided_assembly_ctr_cko_merged_novel_exon_filtered.fa \
-g /mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa \
/mnt/cbis/home/yongshan/at_pipeline/temp_results/stringtie_filtering/exon_level/nestin_ctx_e14_ref_guided_assembly_ctr_cko_merged_novel_exon_filtered.gtf
