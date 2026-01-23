import os
rule sleuth:
	input:
		t2g_aug_collapsed=os.path.join(AT_DIR,"t2g_augment_collapsed.csv"),
		t2g_aug_uncollapsed=os.path.join(AT_DIR,"t2g_augment_uncollapsed.csv"),
		kallisto_tsv=expand(os.path.join(KQ_DIR,"{sample_name}/abundance.tsv"), sample_name = SAMPLES)
	output:
		result_collapsed_tpm=os.path.join(SL_DIR,"collapsed_differential_transcript_analysis_tpm.csv"),
		result_collapsed=os.path.join(SL_DIR,"collapsed_differential_transcript_analysis.csv"),
		result_uncollapsed_tpm=os.path.join(SL_DIR,"uncollapsed_differential_transcript_analysis_tpm.csv"),
		result_uncollapsed=os.path.join(SL_DIR,"uncollapsed_differential_transcript_analysis.csv")
	params:
		kallisto_quant_out=KQ_DIR,
		design_matrix= samples_full_path
	log:
		"logs/sleuth.log"
	conda:
		"../../../../envs/sleuth.yaml"
	threads:
		12
	script:
		"../scripts/sleuth.R"