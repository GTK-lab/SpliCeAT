# Snakemake Parameters
gene_events      <- snakemake@input[["consensus_events"]]
consensus_events <- snakemake@output[["consensus_events_LSVs"]]
consensus_counts <- snakemake@output[["consensus_events_counts"]]
log_file         <- snakemake@log[[1]]

# Logger Setup
library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")
lgr$info("Starting Consensus Analysis...")

lgr$info("Loading libraries...")
# Libraries
suppressMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(tidyr)
  library(gtools)
})

# Helper Function
clean_paste <- function(x) {
  x <- x[!is.na(x) & x != "" & x != "NA"]
  if(length(x) == 0) return(NA_character_)
  paste(unique(x), collapse = ",")
}

# 1. Loading and Indexing
lgr$info("Loading input file...")
all_lsvs <- fread(gene_events) %>%
  as.data.frame() %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end),
         dpsi = as.numeric(dpsi),
         row_idx = row_number())

# Split by Type
junctions <- all_lsvs %>% filter(feature_type == "splice_junction")
exons     <- all_lsvs %>% filter(feature_type == "exon_node")
ir_events <- all_lsvs %>% filter(feature_type == "intron_retention")

# GRanges Objects
j_gr  <- GRanges(junctions$chr, IRanges(junctions$start, junctions$end), junctions$strand)
e_gr  <- GRanges(exons$chr, IRanges(exons$start, exons$end), exons$strand)
ir_gr <- GRanges(ir_events$chr, IRanges(ir_events$start, ir_events$end), ir_events$strand)

# 2. MATCHING LOGIC (maxgap=1)
# 2.1. Range Matching -> J-vs-J and IR-vs-IR
lgr$info("Direct Range matching...")
hits_j  <- as.data.frame(findOverlaps(j_gr, j_gr, type="equal", maxgap=1)) %>% filter(queryHits != subjectHits)
hits_ir <- as.data.frame(findOverlaps(ir_gr, ir_gr, type="equal", maxgap=1)) %>% filter(queryHits != subjectHits)

# 2.2. Boundary Handshakes -> for Junctions flanked by Exons
lgr$info("Boundary Matching...")
j_starts <- GRanges(seqnames(j_gr), IRanges(start(j_gr), width=1), strand(j_gr))
j_ends   <- GRanges(seqnames(j_gr), IRanges(end(j_gr), width=1), strand(j_gr))
e_starts <- GRanges(seqnames(e_gr), IRanges(start(e_gr), width=1), strand(e_gr))
e_ends   <- GRanges(seqnames(e_gr), IRanges(end(e_gr), width=1), strand(e_gr))

h_start <- as.data.frame(findOverlaps(j_starts, e_ends,   maxgap=1))
h_end   <- as.data.frame(findOverlaps(j_ends,   e_starts, maxgap=1))

strict_junctions <- h_start %>%
  mutate(exon_tool = exons$tool[subjectHits]) %>%
  dplyr::inner_join(h_end %>% mutate(exon_tool = exons$tool[subjectHits]),
  			by = c("queryHits", "exon_tool"), relationship = "many-to-many")

# For Exons: Join on Exon ID (subjectHits) AND the Junction Tool (tool)
strict_exons <- h_end %>%
  mutate(junction_tool = junctions$tool[queryHits]) %>%
  dplyr::inner_join(h_start %>% mutate(junction_tool = junctions$tool[queryHits]),
    			by = c("subjectHits", "junction_tool"), relationship = "many-to-many"
  )

# 2.3. IR vs Leafcutter Splice Junction
lgr$info("Matching IR events against Leafcutter Splice Junctions...")
leaf_j_idx <- which(junctions$tool == "Leafcutter")
if(length(leaf_j_idx) > 0){
  ir_leaf_val <- as.data.frame(findOverlaps(ir_gr, j_gr[leaf_j_idx], type="equal", maxgap=1)) %>%
    mutate(ir_row = ir_events$row_idx[queryHits],
           lc_row = junctions$row_idx[leaf_j_idx[subjectHits]],
           ir_dpsi = ir_events$dpsi[queryHits],
           lc_dpsi = junctions$dpsi[leaf_j_idx[subjectHits]],
           ir_tool = ir_events$tool[queryHits]) %>%
	filter(sign(ir_dpsi) != sign(lc_dpsi))
} else { ir_leaf_val <- data.frame() }

# 3. EVIDENCE AGGREGATION
lgr$info("Aggregating evidence across tools...")
evidence_list <- list()
if(nrow(hits_j) > 0)  evidence_list[[1]] <- data.frame(idx = junctions$row_idx[hits_j$queryHits], other_tool = junctions$tool[hits_j$subjectHits])
if(nrow(hits_ir) > 0) evidence_list[[2]] <- data.frame(idx = ir_events$row_idx[hits_ir$queryHits], other_tool = ir_events$tool[hits_ir$subjectHits])
if(nrow(strict_junctions) > 0) evidence_list[[3]] <- data.frame(idx = junctions$row_idx[strict_junctions$queryHits], other_tool = strict_junctions$exon_tool)
if(nrow(strict_exons) > 0)     evidence_list[[4]] <- data.frame(idx = exons$row_idx[strict_exons$subjectHits], other_tool = strict_exons$junction_tool)
if(nrow(ir_leaf_val) > 0) {
  evidence_list[[5]] <- data.frame(idx = ir_leaf_val$ir_row, other_tool = "Leafcutter")
  evidence_list[[6]] <- data.frame(idx = ir_leaf_val$lc_row, other_tool = ir_leaf_val$ir_tool)
}

valid_evidence <- bind_rows(evidence_list) %>% filter(!is.na(idx)) %>% distinct()
validated_rows <- valid_evidence %>% group_by(idx) %>% summarise(n_tools = n_distinct(other_tool) + 1) %>% filter(n_tools >= 2)

lgr$info(sprintf("Total 3-tool agreement events: %d", sum(validated_rows$n_tools == 3)))

# 4. FINAL MATRIX
lgr$info("Formatting final output matrix...")
consensus_matrix <- all_lsvs %>%
  filter(row_idx %in% validated_rows$idx) %>%
  left_join(valid_evidence, by = c("row_idx" = "idx"), relationship = "many-to-many") %>%
  group_by(chr, strand, start, end) %>%
  summarise(
    gene_ids = clean_paste(gene_id),
    gene_names = clean_paste(gene_name),

    # --- MAJIQ ---
    majiq_feature_type = case_when(
        "Majiq" %in% tool ~ feature_type[tool == "Majiq"][1],
        "Majiq" %in% other_tool ~ "(boundary_match)",
        TRUE ~ NA_character_
    )[1],
    majiq_dpsi = clean_paste(dpsi[tool == "Majiq"]),
	majiq_prob = clean_paste(prob[tool == "Majiq"]),

    # --- WHIPPET ---
    whippet_feature_type = case_when(
        "Whippet" %in% tool ~ feature_type[tool == "Whippet"][1],
        "Whippet" %in% other_tool ~ "(boundary_match)",
        TRUE ~ NA_character_
    )[1],
    whippet_dpsi = clean_paste(dpsi[tool == "Whippet"]),
	whippet_prob = clean_paste(prob[tool == "Whippet"]),

    # --- LEAFCUTTER ---
    leafcutter_feature_type = case_when(
      "Leafcutter" %in% tool ~ feature_type[tool == "Leafcutter"][1],
      "Leafcutter" %in% other_tool & any(feature_type == "intron_retention") ~ "(IR_sign_match)",
      "Leafcutter" %in% other_tool ~ "(boundary_match)",
      TRUE ~ NA_character_
    )[1],
    leafcutter_dpsi = clean_paste(dpsi[tool == "Leafcutter"]),
	leafcutter_p_adj = clean_paste(p_adj[tool == "Leafcutter"]),

    n_tools = n_distinct(c(tool, other_tool), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(n_tools >= 2) %>%
  select(chr, strand, start, end, gene_ids, gene_names, n_tools,
         starts_with("majiq"), starts_with("whippet"), starts_with("leafcutter"))
#   arrange(chr, start)

chr_levels <- c(paste0(1:22), "X", "Y", "M", "MT",
                paste0("chr", 1:22), "chrX", "chrY", "chrM", "chrMT")

# 2. Perform the sort using an explicit factor match
consensus_matrix <- consensus_matrix %>%
  mutate(temp_rank = match(as.character(chr), chr_levels)) %>%
  arrange(temp_rank, start) %>%
  select(-temp_rank)

write.table(consensus_matrix, consensus_events, sep="\t", row.names=FALSE, quote=FALSE)

# Generate Gene Summary
gene_summary <- all_lsvs %>%
  filter(row_idx %in% validated_rows$idx) %>%
  group_by(gene_id, gene_name) %>%
  summarise(n_tools = n_distinct(tool),
            total_events = n(),
            .groups = "drop") %>%
  filter(n_tools >= 2) %>%
  arrange(desc(n_tools), desc(total_events))

write.table(gene_summary, consensus_counts, sep="\t", row.names=FALSE, quote=FALSE)
lgr$info(sprintf("Complete. %d consensus events identified.", nrow(consensus_matrix)))