
# Snakemake Parameters
output_consensus <- snakemake@output[["gene_consensus"]]
event_out	     <- snakemake@output[["consensus_events"]]
log_file         <- snakemake@log[[1]]

library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")

# Libraries
lgr$info("Loading libraries...")
suppressMessages({
  library(stringr)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

lgr$info("Loading masterlists...")
majiq_df      <- fread(snakemake@input[["majiq_ms"]]) %>% as.data.frame()
whippet_df    <- fread(snakemake@input[["whippet_ms"]]) %>% as.data.frame()
leafcutter_df <- fread(snakemake@input[["leafcutter_ms"]]) %>% as.data.frame()

lgr$info("Aggregating gene-level detections across tools...")

all_events <- bind_rows(majiq_df, whippet_df, leafcutter_df) %>%
  separate_rows(gene_id, sep = ",") %>%
  mutate(gene_id = str_trim(gene_id)) %>%
  group_by(tool, gene_id, lsv_id, strand) %>%
  slice_max(order_by = abs(dpsi), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  distinct()

gene_consensus <- all_events %>%
	filter(!is.na(gene_id)) %>%
	group_by(gene_id) %>%
	summarise(
		gene_name = paste(unique(na.omit(gene_name)), collapse = ", "),
		n_tools = n_distinct(tool),
		tools_found = paste(sort(unique(tool)), collapse = ", "),
		total_lsv_count = n(),
		.groups = "drop"
	) %>%
	filter(n_tools >= 2) %>%
	arrange(desc(n_tools), gene_id)

lgr$info("Retrieving all LSVs from consensus genes...")
consensus_event_list <- all_events %>%
  filter(gene_id %in% gene_consensus$gene_id) %>%
  mutate(start = as.numeric(start)) %>%
  arrange(chr, start)

consensus_event_list <- consensus_event_list[gtools::mixedorder(consensus_event_list$chr), ]

lgr$info(sprintf("Consensus Analysis Complete."))
lgr$info(sprintf("Genes found by all 3 tools: %d", sum(gene_consensus$n_tools == 3)))
lgr$info(sprintf("Genes found by exactly 2 tools: %d", sum(gene_consensus$n_tools == 2)))

# 5. Save the output
write.table(gene_consensus, output_consensus, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(consensus_event_list, event_out, sep = "\t", row.names = FALSE, quote = FALSE)

lgr$info(sprintf("Consensus Gene stats saved at: %s", output_consensus))
lgr$info(sprintf("Corresponding Gene Events saved at: %s", event_out))