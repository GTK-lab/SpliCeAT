
# Snakemake Parameters
majiq_in       <- snakemake@input[["majiq_in"]]
whippet_in     <- snakemake@input[["whippet_in"]]
leafcutter_in  <- snakemake@input[["leafcutter_in"]]

species        <- snakemake@params[["species"]][[1]] #e.g. Mus_musculus
version        <- snakemake@params[["ensembl"]][[1]]

masterlist_dir <- snakemake@params[["masterlist_dir"]]

if (!dir.exists(masterlist_dir)) dir.create(masterlist_dir, recursive = TRUE)
whippet_out     <- snakemake@output[["whippet_out"]]
majiq_out       <- snakemake@output[["majiq_out"]]
leafcutter_out  <- snakemake@output[["leafcutter_out"]]

# Logger Setup
log_file  <- snakemake@log[[1]]
library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")

# Libraries
lgr$info("Loading libraries...")
suppressMessages({
  library(biomaRt)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# LEAFCUTTER: only splice junctions
lgr$info("LEAFCUTTER: Processing LSVs...")
leafcutter <- read.table(leafcutter_in, sep = "\t", header = TRUE) %>%
  filter(p.adjust < 0.05 & abs(deltapsi) >= 0.2) %>%
  mutate(clean_id = sub("\\.[0-9]+$", "", str_remove(genes, "^gene:")),
         clean_id = ifelse(clean_id == "" | genes == ".", NA_character_, clean_id),
		 start = as.numeric(start),
         end = as.numeric(end),
         tool = "Leafcutter",
		 feature_type = "splice_junction",
         lsv_id = paste0(chr, ":", start, "-", end, ":", strand),
         p_adj = as.numeric(p.adjust)) %>%
  dplyr::select(gene_id=clean_id, chr, strand, start, end, tool, feature_type, lsv_id, dpsi = deltapsi, p_adj)

lgr$info(sprintf("LEAFCUTTER: %d entries are missing gene annotations (NA).", sum(is.na(leafcutter$gene_id))))
lgr$info(sprintf("LEAFCUTTER: Filtering Complete. %d entries present.", nrow(leafcutter)))


# WHIPPET: exons and intron retentions
lgr$info("WHIPPET: Processing LSVs...")
whippet_raw <- readLines(whippet_in)
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
         lsv_id = paste0(chr, ":", start, "-", end, ":", Strand)) %>%
  dplyr::select(gene_id= genes, chr, strand = Strand, start, end, tool, feature_type, lsv_id, dpsi = DeltaPsi, prob = Probability)

lgr$info(sprintf("WHIPPET: Filtering Complete. %d entries present.", nrow(whippet)))


# MAJIQ: only splice junctions
lgr$info("MAJIQ: Processing LSVs....")
majiq_raw <- read.table(majiq_in, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(majiq_raw) <- tolower(colnames(majiq_raw))
majiq <- majiq_raw %>%
  mutate(across(c(mean_dpsi_per_lsv_junction, probability_changing), as.numeric)) %>%
  filter(abs(mean_dpsi_per_lsv_junction) >= 0.2 & probability_changing >= 0.95) %>%
  separate(junctions_coords, into = c("s_raw", "e_raw"), sep = "-", extra = "drop") %>%
  mutate(start = as.integer(pmin(as.numeric(s_raw), as.numeric(e_raw))),
         end = as.integer(pmax(as.numeric(s_raw), as.numeric(e_raw))),
		 tool = "Majiq",
		 feature_type = ifelse(lsv_type_rest == "i", "intron_retention", "splice_junction"),
         lsv_id = paste0(seqid, ":", start, "-", end, ":", strand)) %>%
  group_by(gene_id, seqid, strand, start, end, lsv_id, feature_type) %>%
  slice_max(order_by = abs(mean_dpsi_per_lsv_junction), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(gene_id, chr = seqid, strand, start, end, tool, feature_type, lsv_id,
                dpsi = mean_dpsi_per_lsv_junction, prob = probability_changing) %>%
		mutate(gene_id = sub("\\.[0-9]+$", "", gene_id))

lgr$info(sprintf("MAJIQ: Filtering Complete. %d entries present.", nrow(majiq)))

# ENSEMBL: Gene name fetching

lgr$info("ENSEMBL: Fetch Gene Names...")

all_ids <- unique(na.omit(c(leafcutter$gene_id, whippet$gene_id, majiq$gene_id)))

target_dataset <- switch(
  species,
  "Mus_musculus" = "mmusculus_gene_ensembl",
  "Homo_sapiens" = "hsapiens_gene_ensembl",
  stop(paste0("Unsupported organism: ", species))
)

lgr$info(sprintf("ENSEMBL: Connecting to version %s archive...", version))

ensembl_mart <- tryCatch({
  useEnsembl(biomart = "genes", dataset = target_dataset, version = version)
}, error = function(e) {
  lgr$warn(sprintf("ENSEMBL: Version %s archive unreachable. Switching to live mirror fallback...", version))
  return(NULL)
})

if (is.null(ensembl_mart)) {
  mirrors_to_try <- c("useast", "uswest", "asia")

  for (m in mirrors_to_try) {
    lgr$info(sprintf("ENSEMBL: Connecting to latest live version on mirror: %s", m))
    ensembl_mart <- tryCatch({
      useEnsembl(biomart = "genes", dataset = target_dataset, mirror = m)
    }, error = function(e) { return(NULL) })

    if (!is.null(ensembl_mart)) {
      lgr$info(sprintf("ENSEMBL: Connected to live version via %s", m))
      break
    }
  }
}

if (!is.null(ensembl_mart)) {
  id_to_name <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = all_ids,
    mart       = ensembl_mart
  ) %>% rename(gene_id = ensembl_gene_id, gene_name = external_gene_name)
} else {
  lgr$error("CRITICAL: All Ensembl connection attempts failed.")
  id_to_name <- data.frame(gene_id = character(), gene_name = character())
}

lgr$info("Merging gene names and exporting final masterlists...")

finalize_and_save <- function(df, out_path, tool_name) {
  if (nrow(id_to_name) > 0) {
    df <- df %>% left_join(id_to_name, by = "gene_id")
  } else {
    df$gene_name <- NA
  }

  df <- df %>%
    mutate(start = as.numeric(start), end = as.numeric(end)) %>%
	arrange(chr, start) %>%
    dplyr::select(gene_id, gene_name, everything())

  df <- df[gtools::mixedorder(df$chr), ]

  write.table(df, out_path, sep="\t", row.names=FALSE, quote=FALSE)

  lgr$info(sprintf("MASTERLISTS: %s output list generated at: %s", tool_name, out_path))
}

finalize_and_save(leafcutter, leafcutter_out, "Leafcutter")
finalize_and_save(whippet, whippet_out, "Whippet")
finalize_and_save(majiq, majiq_out, "Majiq")

lgr$info(sprintf("MASTERLISTS generated at %s.", masterlist_dir))

