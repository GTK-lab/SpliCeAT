#transcriptome augmentation
import os

rule prep_for_AT:
	input:
		ds_files = ds_det_output,
		gtf_merged=os.path.join(MG_DIR,"merged_stringtie_assembly.gtf")
	output:
		merged_filtered_gtf=os.path.join(MG_DIR,"merged_stringtie_assembly_novel_exon_filtered.gtf"),
		merged_fil_withRef_gtf=os.path.join(MG_DIR,"merged_stringtie_assembly_novel_exon_filtered_with_reference.gtf")
	log:
		"logs/prep_for_AT.log"
	conda:
		"../../../../envs/for_R.yaml"
	params:
		masterlist_dir = ML_DIR,
		organism = config["ref"]["species"],
		ensembl=config["ref"]["ensembl_release"]
	threads:
		4
	script:
		"../scripts/prep_for_AT.R"

rule collapse_transcripts:
	input:
		gtf_merged=os.path.join(MG_DIR,"merged_stringtie_assembly.gtf"),
		gtf_filtered=os.path.join(MG_DIR,"merged_stringtie_assembly_novel_exon_filtered.gtf")
	output:
		uncollapsed=os.path.join(AT_DIR,"t2g_augment_uncollapsed.csv"),
		collapsed=os.path.join(AT_DIR,"t2g_augment_collapsed.csv")
	log:
		"logs/collapse_transcripts.log"
	conda:
	   "../../../../envs/for_R.yaml"
	params:
		organism = config["ref"]["species"],
		ensembl=config["ref"]["ensembl_release"]
	threads:
		4
	script:
		"../scripts/collapse_transcripts.R"

rule get_novel_sequence:
	input:
		os.path.join(MG_DIR,"merged_stringtie_assembly_novel_exon_filtered.gtf")
	output:
		novel=os.path.join(AT_DIR,"merged_stringtie_assembly_novel_exon_filtered.fa"),
		merged=os.path.join(AT_DIR,"augmented_transcripts.fa")
	params:
		genome=genome_file_path(gz=False,dna=True),
		transcripts=genome_file_path(gz=False,dna=False),
	conda:
		"../../../../envs/genomic_utils.yaml"
	threads:
		4
	shell:
		"gffread -w {output.novel} -g {params.genome} {input} && "
		"cat {params.transcripts} {output.novel} > {output.merged}"