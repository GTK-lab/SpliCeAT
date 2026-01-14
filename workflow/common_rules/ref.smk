## ref.smk
## handle reference filenames and move to common results folder
import os

# config >> resolve pathnames
basepath=config['pipeline_dir'] #base directory full path

sample_file=config["samples"] #design.tsv path wrt config file
samples_full_path = os.path.join(basepath,'config',sample_file)
comparison_groups = config["experiment"]["groups"]

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
	return f"results/{genome_release_name()}{dna_str}.fa{gz_str}"

def gtf_file_path(filtered=False,gz=True):
	gz_str = ".gz" if gz else ""
	filt_str = "_filtered" if filtered else ""
	return f"results/{genome_release_name()}{filt_str}.gtf{gz_str}"

def gff3_file_path(filtered=False,gz=True):
	gz_str = ".gz" if gz else ""
	filt_str = "_filtered" if filtered else ""
	return f"results/{genome_release_name()}{filt_str}.gff3{gz_str}"

def annotation_db_path():
	return f"results/{genome_release_name()}.sqlite3"

# gunzip rule
rule uncompress:
	wildcard_constraints:
		filename=r".+\.(fa|diff|gtf|gff3|gff)"
	input:
		"{filename}.gz"
	conda:
		"../envs/bgzip.yaml"
	output:
		"{filename}",
	wrapper:
		"v5.8.3/bio/bgzip"

