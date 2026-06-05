## ref.smk
## handle reference filenames and move to common results folder
import os

# reference file name
def genome_release_name():
	species=config['ref']['species']
	build=config['ref']['build']
	release=str(config['ref']['ensembl_release'])
	flavor=config['ref']['flavor']
	if flavor:
		release+= "_"
	return f"{species}_{build}_{release}{flavor}"

def genome_file_path(gz=True,dna=True):
	gz_str = ".gz" if gz else ""
	dna_str = "_dna" if dna else "_cds"
	filename = f"{genome_release_name()}{dna_str}.fa{gz_str}"
	return os.path.join(RF_DIR, filename)

def gtf_file_path(filtered=False,gz=True):
	gz_str = ".gz" if gz else ""
	filt_str = "_filtered" if filtered else ""
	filename = f"{genome_release_name()}{filt_str}.gtf{gz_str}"
	return os.path.join(RF_DIR, filename)

def gff3_file_path(filtered=False,gz=True):
	gz_str = ".gz" if gz else ""
	filt_str = "_filtered" if filtered else ""
	filename = f"{genome_release_name()}{filt_str}.gff3{gz_str}"
	return os.path.join(RF_DIR, filename)

def annotation_db_path():
	filename= f"{genome_release_name()}.sqlite3"
	return os.path.join(RF_DIR, filename)

rule uncompress:
	wildcard_constraints:
		filename=r".+\.(fa|diff|gtf|gff3|gff)"
	input:
		"{filename}.gz"
	conda:
		"../envs/bgzip.yaml"
	output:
		"{filename}",
	shell:
		"gunzip -c {input} > {output}"