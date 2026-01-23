import os

rule kallisto_index:
	input:
		os.path.join(AT_DIR,"augmented_transcripts.fa")
	output:
		os.path.join(AT_DIR,"kallisto_index_augmented_transcriptome")
	log:
		"logs/kallisto_index.log"
	conda:
		"../../../../envs/kallisto.yaml"
	threads:
		config["kallisto"]["threads"]
	shell:
		"kallisto index -i {output} {input} > {log} 2>&1"

