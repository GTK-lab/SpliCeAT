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

rule gene_consensus:
	input:
		whippet_ms = os.path.join(ML_DIR, "whippet_lsvs_final.tsv"),
		majiq_ms = os.path.join(ML_DIR, "majiq_lsvs_final.tsv"),
		leafcutter_ms = os.path.join(ML_DIR, "leafcutter_lsvs_final.tsv")
	output:
		gene_consensus = os.path.join(ML_DIR, "consensus_genes_count.tsv"),
		consensus_events = os.path.join(ML_DIR, "consensus_genes_LSV.tsv")
	log:
		"logs/gene_consensus.log"
	conda:
		"../../../../envs/for_R.yaml"
	script:
		"../scripts/gene_consensus.R"

rule granges_consensus:
	input:
		consensus_events = os.path.join(ML_DIR, "consensus_genes_LSV.tsv")
	output:
		consensus_events_counts = os.path.join(ML_DIR, "consensus_granges_count.tsv"),
		consensus_events_LSVs = os.path.join(ML_DIR, "consensus_granges_LSV.tsv")
	log:
		"logs/granges_consensus.log"
	conda:
		"../../../../envs/for_R.yaml"
	script:
		"../scripts/granges_consensus.R"

rule gtf_consensus:
	input:
		consensus_events_LSVs = os.path.join(ML_DIR, "consensus_granges_LSV.tsv"),
		stringtie_gtf=os.path.join(ST_DIR,"stringtie_assembly.gtf"),
		reference_gtf=gtf_file_path(filtered=False,gz=False) #.gtf
	output:
		novel_gtf=os.path.join(AT_DIR,"consensus_novel_transcripts.gtf"),
		augmented_gtf=os.path.join(AT_DIR,"augmented_transcriptome.gtf"),
		validated_lsv_tsv= os.path.join(ML_DIR, "consensus_gtf_LSV.tsv")
	log:
		"logs/gtf_consensus.log"
	conda:
		"../../../../envs/for_R.yaml"
	params:
		species = config["ref"]["species"]
	threads:
		4
	script:
		"../scripts/gtf_consensus.R"

rule collapse_transcripts:
	input:
		gtf_stringtie=os.path.join(ST_DIR,"stringtie_assembly.gtf"),
		gtf_novel_only=os.path.join(AT_DIR,"consensus_novel_transcripts.gtf")
	output:
		uncollapsed=os.path.join(TM_DIR,"t2g_augment_uncollapsed.csv"),
		collapsed=os.path.join(TM_DIR,"t2g_augment_collapsed.csv")
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
		os.path.join(AT_DIR,"consensus_novel_transcripts.gtf")
	output:
		novel=os.path.join(TM_DIR,"consensus_novel_transcripts.fa"),
		merged=os.path.join(TM_DIR,"augmented_transcripts.fa")
	log:
		"logs/get_novel_sequence.log"
	params:
		genome=genome_file_path(gz=False,dna=True),
		transcripts=genome_file_path(gz=False,dna=False),
	conda:
		"../../../../envs/genomic_utils.yaml"
	threads:
		4
	shell:
		"gffread -w {output.novel} -g {params.genome} {input} && "
		"cat {params.transcripts} {output.novel} > {output.merged} 2>> {log}"