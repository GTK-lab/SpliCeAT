# Snakemake Parameters
consensus_events <- snakemake@output[["consensus_events_LSVs"]]
consensus_counts <- snakemake@output[["consensus_events_counts"]]
log_file         <- snakemake@log[[1]]

# Logger Setup
library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")
lgr$info("Starting Consensus Analysis")

# Libraries
lgr$info("Loading libraries...")
suppressMessages({
  library(stringr)
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyr)
  library(lgr)
  library(igraph)
  library(GenomicFeatures)
  library(gtools)
})

# Accessing Snakemake variables
input_events <- fread(snakemake@input[["consensus_events"]]) %>% as.data.frame()

lgr$info("CLUSTERING: Converting outputs to GRanges Objects...")
all_lsvs <- input_events %>%
  mutate(start = as.numeric(start), end = as.numeric(end))

all_gr <- GRanges(all_lsvs$chr, IRanges(all_lsvs$start, all_lsvs$end), all_lsvs$strand)
start_pts <- GRanges(seqnames(all_gr), IRanges(start(all_gr), width=1), strand(all_gr))
end_pts   <- GRanges(seqnames(all_gr), IRanges(end(all_gr), width=1), strand(all_gr))

hits_s <- as.data.frame(findOverlaps(start_pts, start_pts, maxgap = 1))
hits_e <- as.data.frame(findOverlaps(end_pts, end_pts, maxgap = 1))

lgr$info("CLUSTERING: Building Graph of neighbourhood events...")

mirrors <- inner_join(hits_s, hits_e, by = c("queryHits", "subjectHits"))
docking <- rbind(hits_s, hits_e) %>%
  mutate(
    type_q = all_lsvs$feature_type[queryHits],
    type_s = all_lsvs$feature_type[subjectHits]
  ) %>%
  filter(
    (type_q == "exon_node" & type_s == "splice_junction") |
    (type_q == "splice_junction" & type_s == "exon_node")
  )

lgr$info("CLUSTERING: Filtering for repeats...")
valid_connections <- rbind(mirrors[,1:2], docking[,1:2]) %>% unique()
g <- graph_from_edgelist(as.matrix(valid_connections), directed = FALSE)
all_lsvs$cluster <- components(g)$membership

# 3.VALIDATION
lgr$info("VALIDATION: Checking for overlaps between tools...")
site_validation <- rbind(
  data.frame(idx = hits_s$queryHits, tool = all_lsvs$tool[hits_s$subjectHits]),
  data.frame(idx = hits_e$queryHits, tool = all_lsvs$tool[hits_e$subjectHits])
) %>%
  group_by(idx) %>%
  summarise(tools_at_this_site = n_distinct(tool)) %>%
  filter(tools_at_this_site >= 2)

valid_clusters <- all_lsvs %>%
  mutate(row_idx = row_number()) %>%
  filter(row_idx %in% site_validation$idx) %>%
  pull(cluster) %>% unique()

lgr$info("VALIDATION: Generating Final Consensus Matrix...")

clean_paste <- function(x) {
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if(length(x) == 0) return(NA)
  paste(unique(x), collapse = ",")
}

consensus_matrix <- all_lsvs %>%
  filter(cluster %in% valid_clusters) %>%
  group_by(cluster, chr, strand) %>%
  summarise(
	gene_ids = clean_paste(gene_id),
    gene_names = clean_paste(gene_name),

    final_start = min(start),
    final_end = max(end),
    n_tools = n_distinct(tool),

	majiq_type = clean_paste(feature_type[tool == "Majiq"]),
    majiq_dpsi = clean_paste(formatC(dpsi[tool == "Majiq"], digits = 3, format = "f")),
    majiq_prob = clean_paste(formatC(prob[tool == "Majiq"], digits = 3, format = "f")),

    # Tool-specific summaries
    whippet_type = clean_paste(feature_type[tool == "Whippet"]),
    whippet_dpsi = clean_paste(formatC(dpsi[tool == "Whippet"], digits = 3, format = "f")),
    whippet_prob = clean_paste(formatC(prob[tool == "Whippet"], digits = 3, format = "f")),

    leafcutter_event_type = paste(unique(feature_type[tool == "Leafcutter"]), collapse = "|"),
	leaf_dpsi = paste(formatC(dpsi[tool == "Leafcutter"], digits = 3, format = "f"), collapse = "|"),
	leafcutter_p_adj = paste(formatC(p_adj[tool == "Leafcutter"], digits = 3, format = "e"), collapse = "|"),

    .groups = 'drop'
  ) %>%
  filter(n_tools >= 2) %>%
  dplyr::select(
    gene_ids,
	gene_names,
    chr,
    strand,
    start = final_start,
    end = final_end,
    cluster,
    n_tools,
    everything()
  ) %>%
  arrange(chr, start)

consensus_matrix <- consensus_matrix[gtools::mixedorder(consensus_matrix$chr), ]

gene_summary_tsv <- all_lsvs %>%
  filter(cluster %in% valid_clusters) %>%
  separate_rows(gene_id, gene_name, sep = ",") %>%
  mutate(across(c(gene_id, gene_name), str_trim)) %>%
  group_by(gene_id, gene_name) %>%
  summarise(
    n_tools = n_distinct(tool),
    tools_found = paste(sort(unique(tool)), collapse = ", "),
    # This counts how many unique "Splicing Events" (clusters) exist per gene
    total_cluster_count = n_distinct(cluster),
    .groups = "drop"
  ) %>%
  filter(n_tools >= 2) %>%
  arrange(desc(n_tools), gene_id)

write.table(gene_summary_tsv, consensus_counts, sep="\t", row.names=FALSE, quote=FALSE)

write.table(consensus_matrix, consensus_events, sep="\t", row.names=FALSE, quote=FALSE)
lgr$info(sprintf("VALIDATION: Generated Consensus TSV. %d events exported.", nrow(consensus_matrix)))

