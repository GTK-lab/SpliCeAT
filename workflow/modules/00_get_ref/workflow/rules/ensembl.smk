## ensembl.smk
## download relevant reference files from ensembl ftp server

# reference genome fasta
rule get_dna_zipped:
	output:
		genome_file_path(),
	params:
		species=config["ref"]["species"],
		datatype="dna",
		build=config["ref"]["build"],
		release=config["ref"]["ensembl_release"],
	log:
		"logs/get_genome.log",
	cache: "omit-software"
	wrapper:
		"v5.10.0/bio/reference/ensembl-sequence"

# reference transcriptome fasta
rule get_cds_zipped:
	output:
		temp(genome_file_path(dna=False)),
	params:
		species=config["ref"]["species"],
		datatype="cds",
		build=config["ref"]["build"],
		release=config["ref"]["ensembl_release"],
	log:
		"logs/get_transcriptome.log",
	cache: "omit-software"
	wrapper:
		"v5.10.0/bio/reference/ensembl-sequence"

# reference gtf
rule get_gtf_zipped:
	output:
		gtf_file_path(),
	params:
		species=config["ref"]["species"],
		build=config["ref"]["build"],
		release=config["ref"]["ensembl_release"],
		flavor=config["ref"]["flavor"],
	log:
		"logs/get_gtf_gzipped.log",
	cache: "omit-software"
	wrapper:
		"v5.10.0/bio/reference/ensembl-annotation"

# reference gff3
rule get_gff3_zipped:
	output:
		temp(gff3_file_path()),
	params:
		species=config["ref"]["species"],
		build=config["ref"]["build"],
		release=config["ref"]["ensembl_release"],
	log:
		"logs/get_gff3_gzipped.log",
	cache: "omit-software"
	wrapper:
		"v5.10.0/bio/reference/ensembl-annotation"

# indexing gff
rule make_annotation_db:
	input:
		gff3_file_path(),
	output:
		annotation_db_path(),
	log:
		"logs/annotation_db.log",
	conda:
		"../../../../envs/genomic_utils.yaml"
	script:
		"../scripts/make_annotation_db.py"

# filtering downloaded gff
rule make_filtered_gff:
	input:
		annotation_db_path(),
	params:
		biotype=config['experiment']['biotype'],
		max_support_level = config['experiment']['max_transcript_support_level'],
	output:
		gff3_file_path(filtered=True,gz=False)
	conda:
		"../../../../envs/genomic_utils.yaml"
	log:
		"logs/filter_gff.log"
	script:
		"../scripts/filter_annotation_to_gff.py"

# filtered gff -> filtered gff
rule make_filtered_gtf:
	input:
		gff3_file_path(filtered=True,gz=False),
	output:
		gtf_file_path(filtered=True),
	conda:
		"../../../../envs/genomic_utils.yaml",
	log:
		"logs/filter_gtf.log"
	shell:
		"gffread {input} -T | bgzip > {output} 2>{log}"