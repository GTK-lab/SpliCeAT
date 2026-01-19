## move_results.smk
# moving files to common results dir
# common results >> destination file paths

def updated_file(filepath,final_dir):
	base=os.path.basename(filepath)
	return os.path.join(final_dir,base)

def copied_outputs(input_list,final_dir):
	outputs=[]
	for file in input_list:
		outputs.append(updated_file(file,final_dir))
	return outputs

# regarding reference files
REF_FILES = [
	gtf_file_path(filtered=True,gz=False), #.gtf
	gtf_file_path(filtered=True,gz=True), #.gtf.gz
	genome_file_path(gz=True,dna=True), #.dna.fa.gz
	genome_file_path(gz=False,dna=True), #.dna.fa
	genome_file_path(gz=False,dna=False), #cds.fa
	gff3_file_path(filtered=True,gz=False) #.gff3
]
ref_target_dir=os.path.join(PROJECT_ROOT, "results", "references")

rule copy_ref_files:
	input:
		REF_FILES,
	output:
		copied_outputs(REF_FILES,ref_target_dir),
	params:
		ref_dir = ref_target_dir
	log:
		"logs/copy_references.log"
	shell:
		"""
		mkdir -p {params.ref_dir}
		for f in {input}; do
			cp "$f" {params.ref_dir}/
		done
		"""

# regarding ds_detection files
ds_output = [
	f"results/majiq/{'-'.join(GROUPS)}.het.tsv",
	"results/whippet/delta_psi.diff",
	"results/leafcutter/cluster_significance.txt",
	"results/leafcutter/effect_sizes.txt"
]
ds_output_dir=os.path.join(PROJECT_ROOT, "results", "ds_output")

rule DS_GetOutput:
	input:
		ds_output,
	output:
		copied_outputs(ds_output,ds_output_dir),
	params:
		ref_dir = ds_output_dir
	log:
		"logs/ds_output.log"
	shell:
		"""
		mkdir -p {params.ref_dir}
		for f in {input}; do
			cp "$f" {params.ref_dir}/
		done
		"""

# regarding augment_transcriptome files
aug_index = "results/augmented_transcriptome/kallisto_index_augmented_transcriptome"
aug_output = [
	aug_index,
	"results/augmented_transcriptome/t2g_augment_collapsed.csv",
	"results/augmented_transcriptome/t2g_augment_uncollapsed.csv"
]
aug_output_dir=os.path.join(PROJECT_ROOT, "results", "augment_output")

rule Augment_GetOutput:
	input:
		aug_output,
	output:
		copied_outputs(aug_output,aug_output_dir),
	params:
		ref_dir = aug_output_dir
	log:
		"logs/augment_output.log"
	shell:
		"""
		mkdir -p {params.ref_dir}
		for f in {input}; do
			cp "$f" {params.ref_dir}/
		done
		"""
