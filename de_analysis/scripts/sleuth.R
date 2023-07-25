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

lgr$info("Reading config.yaml...")
config <- read_yaml(config_file_path)
lgr$info("Done.")

# configuring h5 files to get rid of trailing transcript info - required if using Gencode annotations.
# If using Ensembl annotations, skip this part
files <- list.files(paste(config$BASE_PATH, "/de_analysis/results/kallisto_quant_out/", sep=""), 
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

############## ANALYSIS 1 - NO COLLAPSING ################

lgr$info("==== Begin Analysis 1: NO collapsing ====")

lgr$info("Configuring experimental design...")

sample_ids <- dir(file.path(paste(config$BASE_PATH, "/de_analysis/results/kallisto_quant_out", sep="")))
lgr$info("Check that your samples are correct:")
sample_id

kal_dirs <- file.path(paste(config$BASE_PATH, "/de_analysis/results/kallisto_quant_out", sep=""), sample_ids, "abundance.h5")
lgr$info("Check that your h5 file paths are correct:")
kal_dirs

lgr$info("Loading in experimental design...")
metadf <- read.table(paste(config$BASE_PATH, "/de_analysis/config/experimental_design.tsv", sep=""), sep = '\t', header = TRUE, colClasses = c("character"))
metadf <- metadf[order(metadf$sample, decreasing=FALSE),]
metadf <- dplyr::mutate(metadf, path = kal_dirs)
lgr$info("Check that the pairing for each sample is correct:")
metadf

metadf$condition <- as.factor(metadf$condition)
metadf$condition <- relevel(metadf$condition, "ctr")

lgr$info("Starting sleuth analysis...")

so <- sleuth_prep(metadf, ~gender+condition, extra_bootstrap_summary = TRUE,
                  gene_mode=FALSE, transformation_function = function(x) log2(x + 0.5)) %>%
      sleuth_fit(~gender, 'reduced') %>%
      sleuth_fit(~gender+condition, 'full') %>%
      sleuth_lrt('reduced', 'full') %>%
      sleuth_wt(which_beta = 'conditiontrtment')

sleuth_tx_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
sleuth_tx_wald <- sleuth_results(so, 'conditiontrtment', 'Wald', show_all = FALSE, pval_aggregate = FALSE)
                  
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

lgr$info("Finished sleuth analysis, saving results...")

write.csv(sleuth_tx_lrt_sleuth_tx_wald_join,paste(config$BASE_PATH, "/de_analysis/results/sleuth/uncollapsed_differential_transcript_analysis.csv", sep=""))

sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
sleuth_matrix <- as.data.frame(sleuth_matrix)
write.csv(sleuth_matrix, paste(config$BASE_PATH, "/de_analysis/results/sleuth/uncollapsed_differential_transcript_analysis_tpm.csv", sep=""))
lgr$info("==== End Analysis 1 ====")

############## ANALYSIS 2 - WITH COLLAPSING ################

lgr$info("==== Begin Analysis 2: WITH collapsing ====")
lgr$info("Starting sleuth analysis...")

so <- sleuth_prep(metadf, ~gender+condition, target_mapping = t2g_augment_collapsed,
                  aggregation_column = 'collapsed_target_id', extra_bootstrap_summary = TRUE,
                  gene_mode=TRUE, transformation_function = function(x) log2(x + 0.5)) %>%
      sleuth_fit(~gender, 'reduced') %>%
      sleuth_fit(~gender+condition, 'full') %>%
      sleuth_lrt('reduced', 'full') %>%
      sleuth_wt(which_beta = 'conditiontrtment')

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
                  
lgr$info("Finished sleuth analysis, saving results...")
write.csv(sleuth_gene_mode_lrt_sleuth_gene_mode_wald_join, paste(config$BASE_PATH, "/de_analysis/results/sleuth/collapsed_differential_transcript_analysis.csv", sep=""))

sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
sleuth_matrix <- as.data.frame(sleuth_matrix)

write.csv(sleuth_matrix, paste(config$BASE_PATH, "/de_analysis/results/sleuth/collapsed_differential_transcript_analysis_tpm.csv", sep=""))
lgr$info("==== End Analysis 2 ====")


