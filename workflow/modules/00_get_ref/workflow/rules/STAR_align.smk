import os
import pandas as pd

rule star_index:
	input:
		fasta=genome_file_path(gz=False,dna=True),
		gtf=gtf_file_path(filtered=False,gz=False),
	output:
		directory(os.path.join(STR_DIR, f"STAR_index_{config['ref']['build']}")),
	message:
		"Generating STAR index"
	threads: config['STAR']['threads'],
	params:
		sjdbOverhang=int(config['experiment']['readlen']) - 1,
		extra=config['STAR']['index_extra'],
	log:
		f"logs/star_index_{config['ref']['build']}.log",
	wrapper:
		"v3.5.3/bio/star/index"

rule star_pe_multi:
	input:
		fq1=get_fq1,
		fq2=get_fq2,
		idx=os.path.join(STR_DIR, f"STAR_index_{config['ref']['build']}"),
	output:
		aln= os.path.join(STR_DIR,"aligned_BAM","{sample}_Aligned.sortedByCoord.out.bam"),
		log="logs/pe/{sample}/Log.out",
		sj="logs/pe/{sample}/SJ.out.tab"
	log:
		"logs/pe/{sample}.log",
	params:
		extra= lambda wildcards: (
			f"--outSAMtype BAM SortedByCoordinate "
			f"--outSAMattrRGline ID:{get_sample_name(wildcards)} "
			f"{config['STAR']['align_extra']} "
		)
	threads: config['STAR']['threads']
	wrapper:
		"v3.5.3/bio/star/align"

rule index_bams:
	input:
		bam=os.path.join(STR_DIR,"aligned_BAM","{sample}_Aligned.sortedByCoord.out.bam"),
	output:
		bai=os.path.join(STR_DIR,"aligned_BAM","{sample}_Aligned.sortedByCoord.out.bam.bai"),
	log:
		"logs/pe/{sample}/index_bam.log",
	threads: 4
	wrapper:
		"v3.5.3/bio/samtools/index"

def STAR_output():
	if config.get("star", {}).get("activate", False):
		return expand(os.path.join(STR_DIR, "aligned_BAM", "{sample}_Aligned.sortedByCoord.out.bam.bai"), sample=SAMPLES)

