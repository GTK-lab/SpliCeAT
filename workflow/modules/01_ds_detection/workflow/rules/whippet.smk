#whippet.smk
#all whippet rules
import os

whippet_final_files = [
	os.path.join(WP_DIR, "merged.bam"),
	os.path.join(WP_DIR, "merged.bam.bai"),
	os.path.join(WP_DIR, "whippet_index.jls"),
	expand(os.path.join(WP_DIR, "quant/{sample}.psi.gz"),sample=SAMPLES),
	os.path.join(WP_DIR, "whippet_delta_psi.diff.gz"),
	os.path.join(WP_DIR,"whippet_delta_psi.diff"),
]

rule whippet:
	input:
		whippet_final_files,
	log:
		"logs/whippet/whippet.log"

rule whippet_merge_bams:
	params:
		bam_groups = config['whippet']['bam_groups'],
	input:
		bam_files = bamfiles_for_groups(config['whippet']['bam_groups']),
	output:
		merged = os.path.join(WP_DIR, "merged.bam"),
		merged_index = os.path.join(WP_DIR, "merged.bam.bai")
	conda:
		"../../../../envs/genomic_utils.yaml"
	log:
		"logs/whippet/merge_bams.log"
	threads: 8,
	shell:
		"""
		samtools merge -o {output.merged} --threads {threads} {input.bam_files} > {log} 2>&1 ;
		samtools index {output.merged} >> {log} 2>&1
		"""

rule whippet_index:
	input:
		genome_fa = genome_file_path(gz=True,dna=True), #dna.fa.gz
		bam =   os.path.join(WP_DIR, "merged.bam"),
		index = os.path.join(WP_DIR, "merged.bam.bai"),
		gtf = gtf_file_path(filtered=True, gz=False) #filtered.gtf.gz
	output:
		os.path.join(WP_DIR, "whippet_index.jls")
	params:
		bam_min_reads = config["whippet"]["bam_min_reads"]
	log:
		"logs/whippet/whippet_index.log"
	conda:
		"../../../../envs/whippet.yaml"
	shell:
		"julia ${{WHIPPET_PATH}}/whippet-index.jl "
		"--fasta {input.genome_fa} "
		"--bam {input.bam} "
		"--gtf {input.gtf} "
		"--bam-min-reads {params.bam_min_reads} "
		"--index {output} > {log} 2>&1;"

rule whippet_quantify:
	input:
		fq1 = get_fq1,
		fq2 = get_fq2,
		index = os.path.join(WP_DIR, "whippet_index.jls"),
	output:
		os.path.join(WP_DIR, "quant/{sample}.psi.gz"),
	params:
		out = lambda w, output: output[0].removesuffix(".psi.gz")
	log:
		"logs/whippet/{sample}_quant.log",
	conda:
		"../../../../envs/whippet.yaml",
	shell:
		"julia ${{WHIPPET_PATH}}/whippet-quant.jl -o {params.out} -x {input.index} {input.fq1} {input.fq2} 2>{log}"

rule whippet_delta:
	input:
		grp1=expand(
			os.path.join(WP_DIR,"quant/{sample}.psi.gz"),
			sample=annot[annot['group'] == GROUPS[0]]['sample_name']),
		grp2=expand(
			os.path.join(WP_DIR,"quant/{sample}.psi.gz"),sample=annot[annot['group'] == GROUPS[1]]['sample_name'])
	params:
		a_list = lambda wildcards, input: ",".join(input.grp1),
		b_list = lambda wildcards, input: ",".join(input.grp2),
		out_prefix = os.path.join(WP_DIR, "whippet_delta_psi")
	output:
		os.path.join(WP_DIR,"whippet_delta_psi.diff.gz"),
	log:
		"logs/whippet/delta_psi.log",
	conda:
		"../../../../envs/whippet.yaml",
	shell:
		"julia ${{WHIPPET_PATH}}/whippet-delta.jl "
		"-a {params.a_list} -b {params.b_list} -o {params.out_prefix} > {log} 2>&1"
