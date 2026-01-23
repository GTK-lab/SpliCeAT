import os
rule stringtie_assembly:
	input:
		get_bam,
		gtf=gtf_file_path(filtered=True,gz=False) #.gtf
	output:
		os.path.join(ST_DIR,"{sample}_ref_guided_assembly.gtf")
	params:
		juncs=config["stringtie"]["min_juncs"],
		strand=config["stringtie"]["strand"]
	conda:
		"../../../../envs/stringtie.yaml"
	log:
		"logs/stringtie/{sample}_stringtie_assembly.log"
	threads:
		config["stringtie"]["threads"]
	shell:
		"""
		stringtie {input} -G {input.gtf} -o {output} -j {params.juncs} -p {threads} --{params.strand} -v 2>&1 | tee {log}.raw | grep --line-buffered -Ei "warning|error|critical|invalid|fail" > {log} || true
		"""
	#verbose mode + filter to catch warnings to primary log file

rule stringtie_merge:
	input:
		expand(os.path.join(ST_DIR,"{sample}_ref_guided_assembly.gtf"), sample=SAMPLES),
		ref_gtf=gtf_file_path(filtered=True,gz=False) #.gtf
	output:
		os.path.join(MG_DIR,"merged_stringtie_assembly.gtf")
	threads:
		config["stringtie"]["threads"]
	conda:
		"../../../../envs/stringtie.yaml"
	log:
		"logs/stringtie/stringtie_merge.log"
	shell:
		"""
		stringtie --merge {input} -G {input.ref_gtf} -i -o {output} -p {threads} -v 2>&1 | tee {log}.raw | grep --line-buffered -Ei "warning|error|critical|invalid|fail" > {log} || true
		"""
	#verbose mode + filter to catch warnings to primary log file
