# Accessing Snakemake variables
species        <- snakemake@params[["organism"]]
version        <- snakemake@params[["ensembl"]]
gtf_merged     <- snakemake@input[["gtf_merged"]]
ds_files       <- snakemake@input[["ds_files"]]
masterlist_dir <- snakemake@params[["masterlist_dir"]]
log_file       <- snakemake@log[[1]]

# Outputs
output_consensus             <- snakemake@output[["consensus_LSVs"]]
merged_filtered_gtf          <- snakemake@output[["merged_filtered_gtf"]]
merged_fil_withRef_gtf       <- snakemake@output[["merged_fil_withRef_gtf"]]

# Masterlist paths
if (!dir.exists(masterlist_dir)) dir.create(masterlist_dir, recursive = TRUE)
whippet_out      <- file.path(masterlist_dir, "whippet_lsvs_final.tsv")
majiq_out        <- file.path(masterlist_dir, "majiq_lsvs_final.tsv")
leafcutter_out   <- file.path(masterlist_dir, "leafcutter_lsvs_final.tsv")
all_lsvs_out     <- file.path(masterlist_dir, "all_lsvs_final.tsv")
output_consensus <- file.path(masterlist_dir,"consensus_LSVs.tsv")

# Logger Setup
library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")
lgr$info("Starting Consensus Analysis")

# Libraries
lgr$info("Loading libraries...")
suppressMessages({
  library(yaml)
  library(biomaRt)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyr)
  library(lgr)
  library(igraph)
})


## 1. TOOL PROCESSING

# LEAFCUTTER: only splice junctions
lgr$info("Leafcutter: Processing LSVs...")
leafcutter <- read.table(ds_files[3], sep = "\t", header = TRUE) %>%
  filter(p.adjust < 0.05 & abs(deltapsi) >= 0.2) %>%
  mutate(genes = str_remove(genes, "^gene:"),
         genes = ifelse(genes == "" | genes == ".", NA_character_, genes),
         tool = "Leafcutter",
		 feature_type = "intron",
         lsv_id = paste0(chr, ":", start, "-", end),
         p_adj = as.numeric(p.adjust)) %>%
  select(genes, chr, strand, start, end, tool, feature_type, lsv_id, dpsi = deltapsi, p_adj)

lgr$info(sprintf("Leafcutter: %d entries are missing gene annotations (NA).", sum(is.na(leafcutter$genes))))
lgr$info(sprintf("Leafcutter: Filtering Complete. %d entries present.", nrow(leafcutter)))

write.table(leafcutter, leafcutter_out, sep="\t", row.names=F, quote=F)
lgr$info(sprintf("Leafcutter: LSV masterlist generated at %s ",leafcutter_out))

# WHIPPET: exons and intron retentions
lgr$info("Whippet: Processing LSVs...")
whippet_raw <- readLines(ds_files[2])
whippet_raw <- whippet_raw[whippet_raw != "" & !grepl("^#", whippet_raw)]
whippet_df <- as.data.frame(do.call(rbind, strsplit(whippet_raw, "\t", fixed = TRUE)))
colnames(whippet_df) <- c("Gene", "Node", "Coord", "Strand", "Type", "Psi_A", "Psi_B", "DeltaPsi", "Probability", "Complexity", "Entropy")

whippet <- whippet_df %>%
  filter(grepl("^-?[0-9.]+", DeltaPsi) & grepl("^[0-9.]+", Probability)) %>%
  mutate(across(c(DeltaPsi, Probability), as.numeric)) %>%
  filter(abs(DeltaPsi) >= 0.2 & Probability >= 0.95 & Type %in% c("CE","AA","AD","RI")) %>%
  separate(Coord, into = c("chr", "coords"), sep = ":") %>%
  separate(coords, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end),
	     genes = str_remove(str_remove(Gene, "^gene:"), "\\.[0-9]+$"),
         genes = ifelse(genes == "" | genes == ".", NA_character_, genes),
         tool = "Whippet",
		 feature_type = ifelse(Type == "RI", "intron_retention", "exon_node"),
         lsv_id = paste0(chr, ":", start, "-", end)) %>%
  select(genes, chr, strand = Strand, start, end, tool, feature_type, lsv_id, dpsi = DeltaPsi, prob = Probability)

lgr$info(sprintf("Whippet: Filtering Complete. %d entries present.", nrow(whippet)))

write.table(whippet, whippet_out, sep="\t", row.names=F, quote=F)
lgr$info(sprintf("Whippet: LSV masterlist generated at %s ",whippet_out))

# MAJIQ: only splice junctions
lgr$info("Majiq: Processing LSVs....")
majiq_raw <- read.table(ds_files[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(majiq_raw) <- tolower(colnames(majiq_raw))
majiq <- majiq_raw %>%
  mutate(across(c(mean_dpsi_per_lsv_junction, probability_changing), as.numeric)) %>%
  filter(abs(mean_dpsi_per_lsv_junction) >= 0.2 & probability_changing >= 0.95) %>%
  separate(junctions_coords, into = c("s_raw", "e_raw"), sep = "-", extra = "drop") %>%
  mutate(start = as.integer(pmin(as.numeric(s_raw), as.numeric(e_raw))),
         end = as.integer(pmax(as.numeric(s_raw), as.numeric(e_raw))),
         genes = str_remove(gene_id, "^gene:"),
		 tool = "Majiq",
		 feature_type = "intron",
         lsv_id = paste0(seqid, ":", start, "-", end)) %>%
  select(genes, chr = seqid, strand, start, end, tool, feature_type, lsv_id, dpsi = mean_dpsi_per_lsv_junction, prob = probability_changing)

lgr$info(sprintf("Majiq: Filtering Complete. %d entries present.", nrow(majiq)))

write.table(majiq, majiq_out, sep="\t", row.names=F, quote=F)
lgr$info(sprintf("Majiq: LSV masterlist generated at %s.",majiq_out))

## 2. CLUSTERING

lgr$info("CLUSTERING: Converting outputs to GRanges Objects...")
all_lsvs <- bind_rows(leafcutter, whippet, majiq) %>%
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
    (type_q == "exon_node" & type_s == "intron") |
    (type_q == "intron" & type_s == "exon_node")
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

consensus_matrix <- all_lsvs %>%
  filter(cluster %in% valid_clusters) %>%
  group_by(cluster, chr, strand) %>%
  summarise(
    genes = paste(unique(na.omit(genes)), collapse = ","),
    final_start = min(start),
    final_end = max(end),

    n_tools = n_distinct(tool),

	majiq_event_type = paste(unique(feature_type[tool == "Majiq"]), collapse = "|"),
	majiq_dpsi = paste(formatC(dpsi[tool == "Majiq"], digits = 3, format = "f"), collapse = "|"),
	majiq_prob = paste(formatC(prob[tool == "Majiq"], digits = 3, format = "f"), collapse = "|"),

    # Tool-specific summaries
    whippet_event_type = paste(unique(feature_type[tool == "Whippet"]), collapse = "|"),
    whippet_dpsi = paste(formatC(dpsi[tool == "Whippet"], digits = 3, format = "f"), collapse = "|"),
	whippet_prob   = paste(formatC(prob[tool == "Whippet"], digits = 3, format = "f"), collapse = "|"),

    leafcutter_event_type = paste(unique(feature_type[tool == "Leafcutter"]), collapse = "|"),
	leaf_dpsi = paste(formatC(dpsi[tool == "Leafcutter"], digits = 3, format = "f"), collapse = "|"),
	leafcutter_p_adj = paste(formatC(p_adj[tool == "Leafcutter"], digits = 3, format = "e"), collapse = "|"),

    .groups = 'drop'
  ) %>%
  filter(n_tools >= 2) %>%
  select(
    genes,
    chr,
    strand,
    start = final_start,
    end = final_end,
    cluster,
    n_tools,
    everything()
  ) %>%
  mutate(across(where(is.character), ~na_if(., "")))

write.table(consensus_matrix, output_consensus, sep="\t", row.names=FALSE, quote=FALSE)
lgr$info(sprintf("VALIDATION: Generated Consensus TSV. %d events exported.", nrow(consensus_matrix)))

# 4. GTF FILTERING
lgr$info("CONSTRUCTING GTF: Processing merged Stringtie GTF file...")
gtf_full <- rtracklayer::import(gtf_merged)

consensus_gr <- makeGRangesFromDataFrame(
  consensus_matrix,
  keep.extra.columns = TRUE,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end"
)

hits <- findOverlaps(gtf_full, consensus_gr)
overlapping_tx_ids <- unique(mcols(gtf_full[queryHits(hits)])$transcript_id)
novel_supported_tx <- overlapping_tx_ids[grepl("MSTRG", overlapping_tx_ids)]

# novel isoforms only -> merged_stringtie_assembly_novel_exon_filtered.gtf
gtf_novel_only <- gtf_full[mcols(gtf_full)$transcript_id %in% novel_supported_tx]

rtracklayer::export(gtf_novel_only, merged_filtered_gtf, format="gtf")
lgr$info(sprintf("CONSTRUCTING GTF: Filtering complete. %d entries present.", length(novel_supported_tx)))
lgr$info(sprintf("CONSTRUCTING GTF: Novel-only GTF generated at %s", merged_filtered_gtf))

# augmented to reference GTF
ref_pattern <- ifelse(species == "Mus_musculus", "ENSMUST", "ENST")

gtf_with_ref <- c(
  gtf_full[grepl(ref_pattern, mcols(gtf_full)$transcript_id)],
  gtf_novel_only
)

rtracklayer::export(gtf_with_ref, merged_fil_withRef_gtf, format="gtf")
lgr$info("CONSTRUCTING GTF: Complete Augmented GTF generated at %s", merged_fil_withRef_gtf)