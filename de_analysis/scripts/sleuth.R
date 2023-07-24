# sleuth for differential transcript expression analysis

#### CHANGE THIS ####
config_file_path <- "/mnt/cbis/home/yongshan/SpliCeAT/de_analysis/config/config.yaml"
####################

library(lgr)

lgr$info("Loading libraries...")
suppressMessages({
library(rhdf5)
})
lgr$info("Done.")

# load in config containing file paths and other params
lgr$info("Reading config.yaml...")
config <- read_yaml(config_file_path)
lgr$info("Done.")

# configuring h5 files to get rid of trailing transcript info - required if using Gencode annotations. If using Ensembl annotations, skip this part
files <- list.files(paste(config$BASE_PATH, "/de_analysis/kallisto_quant_out/", sep=""), 
                    pattern=".h5", 
                    recursive=TRUE, 
                    full.names=TRUE)

for (currentFile in files) {
  oldids <- h5read(currentFile, "/aux/ids")
  newids <- gsub("\\|.*", "", oldids)
  h5write(newids, currentFile, "/aux/ids")
}

lgr$info("Loading t2g dataframes...")
t2g_augment_uncollapsed <- read.csv(paste(config$BASE_PATH,"/augment_transcriptome/results/augmented_transcriptome/t2g_augment_uncollapsed.csv", sep=""))
t2g_augment_collapsed <- read.csv(paste(config$BASE_PATH,"/augment_transcriptome/results/augmented_transcriptome/t2g_augment_collapsed.csv", sep=""))

############## ANALYSIS 2 - WITH COLLAPSING ################
# follow sleuth walkthrough to load the experimental tsv table
metadf <- read.table("/mnt/cbis/home/yongshan/collapse_transcript_test/experimental_design.tsv", sep = '\t', header = TRUE, colClasses = c("character"))

metadf$condition <- as.factor(metadf$condition)
metadf$condition <- relevel(metadf$condition, "ctr")

head(t2g_augment)
tail(t2g_augment)

so <- sleuth_prep(metadf, ~gender+condition, target_mapping = t2g_augment,
                  aggregation_column = 'collapsed_target_id', extra_bootstrap_summary = TRUE,
                  gene_mode=TRUE, transformation_function = function(x) log2(x + 0.5)) ## NOT BETA but LOG2FC

so <- sleuth_fit(so, ~gender, 'reduced')

so <- sleuth_fit(so, ~gender+condition, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')

so <- sleuth_wt(so,which_beta = 'conditiontrtment')

sleuth_gene_mode_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE, gene_mode = TRUE)
sleuth_gene_mode_wald <- sleuth_results(so, 'conditiontrtment', 'Wald', show_all = FALSE, pval_aggregate = FALSE, gene_mode = TRUE)

colnames(sleuth_gene_mode_lrt) <- c("ens_gene",
                                    "ext_gene",
                                    "target_id",
                                    "pval_lrt",
                                    "qval_lrt",
                                    "test_stat_lrt",
                                    "rss_lrt",
                                    "degrees_free_lrt",
                                    "mean_obs_lrt",
                                    "var_obs_lrt",
                                    "tech_var_lrt",
                                    "sigma_sq_lrt",
                                    "smooth_sigma_sq_lrt",
                                    "final_sigma_sq_lrt")

colnames(sleuth_gene_mode_wald) <- c("ens_gene",
                                     "ext_gene",
                                     "target_id",
                                     "pval_wald",
                                     "qval_wald",
                                     "b_wald",
                                     "se_b_wald",
                                     "mean_obs_wald",
                                     "var_obs_wald",
                                     "tech_var_wald",
                                     "sigma_sq_wald",
                                     "smooth_sigma_sq_wald",
                                     "final_sigma_sq_wald")

sleuth_gene_mode_lrt_sleuth_gene_mode_wald_join <- full_join(sleuth_gene_mode_lrt,sleuth_gene_mode_wald,
                                                             by=c("ens_gene", "ext_gene","target_id"))

write.csv(sleuth_gene_mode_lrt_sleuth_gene_mode_wald_join,"/mnt/cbis/home/yongshan/collapse_transcript_test/tdp43_nestin_ctx_e14/nestin_ctx_e14_gene_mode_lrt_wald_join_augmented_collapsed.csv")

sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
#sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'est_counts')
head(sleuth_matrix) # look at first 5 transcripts, sorted by name
tail(sleuth_matrix)
sleuth_matrix <- as.data.frame(sleuth_matrix)

write.csv(sleuth_matrix,"/mn/cbis/home/yongshan/collapse_transcript_test/tdp43_nestin_ctx_e14/nestin_ctx_e14_tpm_counts_augmented_collapsed.csv")

############## NORMAL DTE ANALYSIS - no collapsing ################
metadf <- read.table("/mnt/cbis/home/yongshan/collapse_transcript_test/experimental_design.tsv", sep = '\t', header = TRUE, colClasses = c("character"))

metadf$condition <- as.factor(metadf$condition)
metadf$condition <- relevel(metadf$condition, "ctr")

so <- sleuth_prep(metadf, ~gender+condition, extra_bootstrap_summary = TRUE,
                  gene_mode=FALSE, transformation_function = function(x) log2(x + 0.5)) ## NOT BETA but LOG2FC

so <- sleuth_fit(so, ~gender, 'reduced')

so <- sleuth_fit(so, ~gender+condition, 'full')

so <- sleuth_lrt(so, 'reduced', 'full')

so <- sleuth_wt(so,which_beta = 'conditiontrtment')

sleuth_tx_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)

sleuth_tx_wald <- sleuth_results(so, 'conditiontrtment', 'Wald', show_all = FALSE, pval_aggregate = FALSE)

#sleuth_live(so)
sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
#sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'est_counts')
head(sleuth_matrix) # look at first 5 transcripts, sorted by name
write.csv(sleuth_matrix,"/mnt/cbis/home/yongshan/collapse_transcript_test/tdp43_nestin_ctx_e14/nestin_ctx_e14_tpm_counts_augmented.csv")


nrow(sleuth_tx_lrt) # 80500
head(sleuth_tx_lrt)
colnames(sleuth_tx_lrt) <- c("target_id",
                             "pval_lrt",
                             "qval_lrt",
                             "test_stat_lrt",
                             "rss_lrt",
                             "degrees_free_lrt",
                             "mean_obs_lrt",
                             "var_obs_lrt",
                             "tech_var_lrt",
                             "sigma_sq_lrt",
                             "smooth_sigma_sq_lrt",
                             "final_sigma_sq_lrt")

nrow(sleuth_tx_wald) # 65224
head(sleuth_tx_wald)
colnames(sleuth_tx_wald) <- c("target_id",
                              "pval_wald",
                              "qval_wald",
                              "b_wald",
                              "se_b_wald",
                              "mean_obs_wald",
                              "var_obs_wald",
                              "tech_var_wald",
                              "sigma_sq_wald",
                              "smooth_sigma_sq_wald",
                              "final_sigma_sq_wald")

sleuth_tx_lrt_sleuth_tx_wald_join <- full_join(sleuth_tx_lrt,sleuth_tx_wald,by=c("target_id"))

head(sleuth_tx_lrt_sleuth_tx_wald_join)
nrow(sleuth_tx_lrt_sleuth_tx_wald_join)

write.csv(sleuth_tx_lrt_sleuth_tx_wald_join,"/mnt/cbis/home/yongshan/collapse_transcript_test/tdp43_nestin_ctx_e14/nestin_ctx_e14_gene_mode_lrt_wald_join_augmented.csv")

