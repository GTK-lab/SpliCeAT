#irfinder.smk
#all irfinder rules
import os

STAR_INDEX_DIR = os.path.join(STR_DIR, f"STAR_index_{config['ref']['build']}")
IRF_STAR_INDEX = STAR_INDEX_DIR if star_active else os.path.join(IR_DIR, "star_index")

IRF_REF_DIR = os.path.join(IR_DIR, "ref", genome_release_name())

IRF_STRANDED = str(config["experiment"]["strandedness"]).lower() in ("yes", "forward", "fr", "rf", "reverse")
IRF_RESULT = "IRFinder-IR-dir.txt" if IRF_STRANDED else "IRFinder-IR-nondir.txt"

IRF_GROUPS = prepare_irfinder_names(GROUPS, sample_for_group)

# STAR required to built IRFinder reference
if not star_active:
	rule irfinder_star_index:
		input:
			fasta=genome_file_path(gz=False, dna=True),
			gtf=gtf_file_path(filtered=False, gz=False),
		output:
			directory(os.path.join(IR_DIR, "star_index")),
		message:
			"Building STAR index for IRFinder-S (STAR alignment is disabled)"
		params:
			overhang=int(config["experiment"]["readlen"]) - 1,
		threads:
			16
		log:
			"logs/irfinder/star_index.log",
		conda:
			"../../../../envs/irfinder.yaml"
		shell:
			"mkdir -p {output} && "
			"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output} "
			"--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} "
			"--sjdbOverhang {params.overhang} --outFileNamePrefix {output}/ > {log} 2>&1"

rule irfinder_build_ref:
	input:
		star_index=IRF_STAR_INDEX,
		fasta=genome_file_path(gz=False, dna=True),
		gtf=gtf_file_path(filtered=False, gz=False), 
	output:
		ref_dir=directory(IRF_REF_DIR),
	message:
		"Generating IRFinder reference"
	params:
		read_length=config["experiment"]["readlen"],
	threads: 
		8
	log:
		"logs/irfinder/build_ref.log",
	conda:
		"../../../../envs/irfinder.yaml"
	shell:
		"IRFinder BuildRefFromSTARRef -l -r {output.ref_dir} -x {input.star_index} "
		"-f {input.fasta} -g {input.gtf} -t {threads} -n {params.read_length} > {log} 2>&1"

rule irfinder_bam:
	input:
		bam=get_bam,
		bai=lambda wildcards: get_bam(wildcards) + ".bai",
		ref=IRF_REF_DIR,
	output:
		result=os.path.join(IR_DIR, "quant", "{sample}", IRF_RESULT),
	params:
		outdir=lambda wildcards, output: os.path.dirname(output.result),
		cnn="" if config["irfinder"]["cnn_validation"] else "-w 0", # -w 0 disables the CNN validation step,
	threads: 
		1 
	log:
		"logs/irfinder/{sample}_bam.log",
	conda:
		"../../../../envs/irfinder.yaml"
	shell:
		"rm -rf {params.outdir};" # clean previous output to allow sample reruns
		"IRFinder BAM -r {input.ref} -d {params.outdir} {params.cnn} {input.bam} > {log} 2>&1"

rule irfinder_diff:
	input:
		grp1=expand(os.path.join(IR_DIR, "quant", "{sample}", IRF_RESULT), sample=sample_for_group(GROUPS[0])),
		grp2=expand(os.path.join(IR_DIR, "quant", "{sample}", IRF_RESULT), sample=sample_for_group(GROUPS[1])),
	output:
		all_results=os.path.join(IR_DIR, "diff", "all_results_DESeq2.tsv"),
	params:
		workdir=IR_DIR,
		grp1=lambda _: f"-g:{IRF_GROUPS[0]} {irfinder_group_files(GROUPS[0])}",
		grp2=lambda _: f"-g:{IRF_GROUPS[1]} {irfinder_group_files(GROUPS[1])}",
		min_irratio=config["irfinder"]["min_irratio"],
		warning_level=config["irfinder"]["warning_level"],
	log:
		"logs/irfinder/diff.log",
	conda:
		"../../../../envs/irfinder.yaml"
	shell:
		"(cd {params.workdir} && IRFinder Diff {params.grp1} {params.grp2} "
		"-o diff -m deseq -ir {params.min_irratio} -wl {params.warning_level}) > {log} 2>&1"

rule irfinder_lsvs:
	input:
		all_results=os.path.join(IR_DIR, "diff", "all_results_DESeq2.tsv"),
		ir_result=os.path.join(IR_DIR, "quant", sample_for_group(GROUPS[0])[0], IRF_RESULT), 
	output:
		lsvs=os.path.join(IR_DIR, "irfinder_lsvs.tsv"),
	params:
		groups=GROUPS[:2],
		irf_groups=IRF_GROUPS,
		samples=sample_for_group(GROUPS[0]) + sample_for_group(GROUPS[1]),
	log:
		"logs/irfinder/lsvs.log",
	conda:
		"../../../../envs/base.yaml"
	script:
		"../scripts/irfinder_lsvs.py"
