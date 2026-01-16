rule stringtie_assembly:
	input:
		get_bam
	output:
		"results/stringtie_assemblies/{sample}_ref_guided_assembly.gtf"
	params:
		gtf=updated_file(gtf_file_path(filtered=True,gz=False),ref_target_dir), #.gtf
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
		stringtie {input} -G {params.gtf} -o {output} -j {params.juncs} -p {threads} --{params.strand} -v 2>&1 | tee {log}.raw | grep --line-buffered -Ei "warning|error|critical|invalid|fail" > {log} || true
		"""
	#verbose mode + filter to catch warnings to primary log file

rule stringtie_merge:
	input:
		expand("results/stringtie_assemblies/{sample}_ref_guided_assembly.gtf", sample=SAMPLES)
	output:
		"results/merged_assembly/merged_stringtie_assembly.gtf"
	params:
		gtf=updated_file(gtf_file_path(filtered=True,gz=False),ref_target_dir) #.gtf
	threads:
		config["stringtie"]["threads"]
	conda:
		"../../../../envs/stringtie.yaml"
	log:
		"logs/stringtie/stringtie_merge.log"
	shell:
		"""
		stringtie --merge {input} -G {params.gtf} -i -o {output} -p {threads} -v 2>&1 | tee {log}.raw | grep --line-buffered -Ei "warning|error|critical|invalid|fail" > {log} || true
		"""
	#verbose mode + filter to catch warnings to primary log file
