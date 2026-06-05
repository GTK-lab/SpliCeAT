#update!!
import os

def kallisto_strandedness(strandedness=config["experiment"]["strandedness"]):
	val = str(strandedness).lower()
	lookup = {
		"yes": "--fr-stranded",
		"forward": "--fr-stranded",
		"fr": "--fr-stranded",
		"rf": "--rf-stranded",
		"reverse": "--rf-stranded",
		"none" : "",
		"no" : "",
		"unstranded": "",
	}
	return lookup.get(val, "")

rule kallisto_test:
	input:
		kallisto_tsv=expand(os.path.join(KQ_DIR,"{sample_name}/abundance.tsv"), sample_name = SAMPLES)

rule kallisto_quant:
	input:
		unpack(get_fastq_files)
	output:
		os.path.join(KQ_DIR,"{sample}/abundance.tsv")
	params:
		INDEX =  os.path.join(TM_DIR,"kallisto_index_augmented_transcriptome"),
		OUT_FILE = os.path.join(KQ_DIR,"{sample}"),
		STRAND = kallisto_strandedness(),
		BOOTSTRAPS = config["kallisto"]["bootstraps"],
	threads:
		config["kallisto"]["threads"]
	conda:
		"../../../../envs/kallisto.yaml"
	log:
		"logs/kallisto_quant/{sample}.log"
	message:
		"--- kallisto quant ---"
	shell:
		"kallisto quant -i {params.INDEX} -o {params.OUT_FILE} {params.STRAND} --bootstrap-samples={params.BOOTSTRAPS} --threads={threads} {input} > {log} 2>&1"
