#transcriptome augmentation
import os

rule generate_masterlists:
	input:
		majiq_in = os.path.join(MJ_DIR,f"majiq_expanded_{'-'.join(GROUPS)}.deltapsi.tsv"),
		whippet_in = os.path.join(WP_DIR, "whippet_delta_psi.diff"),
		leafcutter_in = os.path.join(LC_DIR, "leafcutter_lsvs.tsv")
	output:
		whippet_out = os.path.join(ML_DIR, "whippet_lsvs_final.tsv"),
		majiq_out = os.path.join(ML_DIR, "majiq_lsvs_final.tsv"),
		leafcutter_out = os.path.join(ML_DIR, "leafcutter_lsvs_final.tsv")
	params:
		masterlist_dir = ML_DIR,
		species = config["ref"]["species"],
		ensembl = config["ref"]["ensembl_release"]
	log:
		"logs/generate_masterlists.log"
	conda:
		"../../../../envs/for_R.yaml"
	script:
		"../scripts/generate_masterlists.R"

rule get_gene_overlap:
	input:
		whippet_ms = os.path.join(ML_DIR, "whippet_lsvs_final.tsv"),
		majiq_ms = os.path.join(ML_DIR, "majiq_lsvs_final.tsv"),
		leafcutter_ms = os.path.join(ML_DIR, "leafcutter_lsvs_final.tsv")
	output:
		gene_consensus = os.path.join(ML_DIR, "Consensus_Genes_count.tsv"),
		consensus_events = os.path.join(ML_DIR, "Consensus_Gene_LSVs.tsv")
	log:
		"logs/gene_consensus.log"
	conda:
		"../../../../envs/for_R.yaml"
	script:
		"../scripts/gene_consensus.R"

rule get_granges_overlap:
	input:
		consensus_events = os.path.join(ML_DIR, "Consensus_Gene_LSVs.tsv")
	output:
		consensus_events_counts = os.path.join(ML_DIR, "Consensus_GRanges_count.tsv"),
		consensus_events_LSVs = os.path.join(ML_DIR, "Consensus_GRanges_LSVs.tsv")
	log:
		"logs/granges_consensus.log"
	conda:
		"../../../../envs/for_R.yaml"
	script:
		"../scripts/LSV_consensus.R"

rule prep_for_AT:
	input:
		consensus_events_LSVs = os.path.join(ML_DIR, "Consensus_GRanges_LSVs.tsv"),
		stringtie_gtf=os.path.join(MG_DIR,"merged_stringtie_assembly.gtf")
	output:
		merged_filtered_gtf=os.path.join(MG_DIR,"stringtie_filtered_novel_exons.gtf"),
		merged_fil_withRef_gtf=os.path.join(MG_DIR,"stringtie_filtered_with_reference.gtf")
	log:
		"logs/prep_for_AT.log"
	conda:
		"../../../../envs/for_R.yaml"
	params:
		species = config["ref"]["species"]
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