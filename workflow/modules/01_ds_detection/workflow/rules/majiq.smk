#majiq.smk
#all majiq rules
import os

majiq_sj_files = [os.path.join(MJ_DIR,f"{f}.sj") for f in annot['bam_stem'] ]
majiq_majiq_files = [os.path.join(MJ_DIR,f"{f}.majiq") for f in annot['bam_stem'] ]

group_prefix = "-".join(GROUPS)
majiq_final_files = [
	expand(os.path.join(MJ_DIR,"{sample}.sj"), sample=annot['bam_stem']),
	os.path.join(MJ_DIR,f"majiq_{group_prefix}.deltapsi.tsv"),
	os.path.join(MJ_DIR,f"majiq_{group_prefix}.het.tsv"),
	os.path.join(MJ_DIR,f"majiq_expanded_{group_prefix}.deltapsi.tsv"),
	os.path.join(MJ_DIR,f"majiq_{group_prefix}.deltapsi.voila"),
	expand(os.path.join(MJ_DIR,"{sample}.majiq"), sample=annot['bam_stem'])
]

rule majiq:
	input:
		majiq_final_files
	output:
		os.path.join(MJ_DIR,".donestamp")
	log:
		"logs/majiq.log"
	conda:
		"../../../../envs/base.yaml"
	shell:
		"touch {output}"

rule majiq_conf:
	input:
		design_tsv=samples_full_path, #full path of design.tsv
		bam_files=bam_files,
	log:
		"logs/majiq/build_ini.log",
	output:
		touch(os.path.join(MJ_DIR,"majiq.ini")),
	conda:
		"../../../../envs/base.yaml"
	script:
		"../scripts/build_majiq_ini.py"

rule majiq_build:
	input:
		gff3=gff3_file_path(filtered=True,gz=False), #.gff3 ##update!!
		conf_file = os.path.join(MJ_DIR,"majiq.ini"),
		bams = bam_files,
	params:
		license=config["majiq"]["license"],
		extra=config["majiq"]["build_extra"],
		output_dir= MJ_DIR
	output:
		os.path.join(MJ_DIR,"splicegraph.sql"),
		temp(majiq_sj_files),
		temp(majiq_majiq_files)
	log:
		"logs/majiq/majiq_build.log",
	conda:
		"../../../../envs/majiq.yaml"
	threads:
		8
	shell:
		"majiq --license {params.license} build {input.gff3} "
		"--logger {log} {params.extra} "
		"-c {input.conf_file} -j {threads} -o {params.output_dir} >/dev/null 2>&1"

rule majiq_delta_psi:
	input:
		sj_files=majiq_sj_files,
		majiq_files=majiq_majiq_files,
		splicegraph=os.path.join(MJ_DIR,"splicegraph.sql"),
	output:
		tsv=os.path.join(MJ_DIR,f"majiq_{group_prefix}.deltapsi.tsv"),
		voila=os.path.join(MJ_DIR,f"majiq_{group_prefix}.deltapsi.voila"),
	log:
		"logs/majiq/majiq_deltapsi.log",
	params:
		output_dir=MJ_DIR,
		license=config['majiq']['license'] ,
		names=lambda _: " ".join(GROUPS),
		raw_prefix=os.path.join(MJ_DIR, group_prefix),
		grp1= lambda _: " ".join(majiq_files(GROUPS[0])),
		grp2= lambda _: " ".join(majiq_files(GROUPS[1])),
	threads:
		16
	conda:
		"../../../../envs/majiq.yaml"
	shell:
		"majiq --license {params.license} deltapsi "
		"-grp1 {params.grp1} -grp2 {params.grp2} "
		"--logger {log} "
		"-j {threads} -o {params.output_dir} -n {params.names} > /dev/null 2> {log}.err ;"
		"mv {params.raw_prefix}.deltapsi.tsv {output.tsv};"
		"mv {params.raw_prefix}.deltapsi.voila {output.voila}"

rule majiq_heterogen:
	input:
		sj_files=majiq_sj_files,
		majiq_files=majiq_majiq_files,
		splicegraph=os.path.join(MJ_DIR,"splicegraph.sql"),
	output:
		voila= os.path.join(MJ_DIR, f"majiq_{group_prefix}.het.voila"),
	log:
		"logs/majiq/majiq_heterogen.log",
	params:
		output_dir=MJ_DIR,
		license=config['majiq']['license'] ,
		names=lambda _: " ".join(GROUPS),
		raw_prefix=os.path.join(MJ_DIR, group_prefix),
		grp1= lambda _: " ".join(majiq_files(GROUPS[0])),
		grp2= lambda _: " ".join(majiq_files(GROUPS[1])),
	threads:
		16
	conda:
		"../../../../envs/majiq.yaml"
	shell:
		"majiq --license {params.license} heterogen "
		"-grp1 {params.grp1} -grp2 {params.grp2} "
		"--logger {log} "
		"-j {threads} -o {params.output_dir} -n {params.names} > /dev/null 2> {log}.err;"
		"mv {params.raw_prefix}.het.voila {output.voila}"

rule majiq_voila_to_tsv:
	input:
		voila=os.path.join(MJ_DIR,f"majiq_{group_prefix}.het.voila"),
		splicegraph=os.path.join(MJ_DIR,"splicegraph.sql"),
	output:
		os.path.join(MJ_DIR,f"majiq_{group_prefix}.het.tsv"),
	log:
		"logs/majiq/het_voila_to_tsv.log",
	conda:
		"../../../../envs/majiq.yaml"
	shell:
		"voila tsv {input.splicegraph} {input.voila} -f {output} > {log} 2>&1"

rule majiq_explode:
	input:
		os.path.join(MJ_DIR,f"majiq_{group_prefix}.deltapsi.tsv"),
	output:
		os.path.join(MJ_DIR,f"majiq_expanded_{group_prefix}.deltapsi.tsv"),
	conda:
		"../../../../envs/base.yaml"
	log:
		"logs/majiq/explode.log"
	script:
		"../scripts/explode_majiq_tsv.py"