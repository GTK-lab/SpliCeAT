# Snakemake Parameters
consensus_events <- snakemake@input[["consensus_events_LSVs"]]
stringtie_gtf    <- snakemake@input[["stringtie_gtf"]]
reference_gtf    <- snakemake@input[["reference_gtf"]]

species                <- snakemake@params[["species"]]
merged_filtered_gtf    <- snakemake@output[["merged_filtered_gtf"]]
merged_fil_withRef_gtf <- snakemake@output[["merged_fil_withRef_gtf"]]
log_file               <- snakemake@log[[1]]

# Logger Setup
library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")

# Libraries
suppressMessages({
    library(rtracklayer)
    library(GenomicFeatures)
    library(txdbmaker)
    library(data.table)
    library(dplyr)
})

# Remove ensembl version numbers and clean metadata columns for matching with LSV file
clean_ids <- function(x) {
    x <- gsub("transcript:|gene:", "", as.character(x))
    x <- ifelse(grepl("^MSTRG", x), x, sub("\\..*$", "", x))
    return(x) }

# 1. TRUTH SET (Ref + Consensus)
lgr$info("TRUTH SET: Extracting introns from Reference GTF...")
ref_txdb <- txdbmaker::makeTxDbFromGFF(reference_gtf)
ref_introns <- GenomicFeatures::intronsByTranscript(ref_txdb) %>% unlist() %>% unique()

lgr$info("TRUTH SET: Extracting junctions from Consensus LSVs...")
consensus_matrix <- fread(consensus_events) %>% as.data.frame()
consensus_ranges <- GRanges(
    seqnames = consensus_matrix$chr,
    ranges = IRanges(start = consensus_matrix$start, end = consensus_matrix$end),
    strand = consensus_matrix$strand
)

safe_ranges <- unique(c(ref_introns, consensus_ranges))

# 2. STRINGTIE GTF
lgr$info("STRINGTIE: Building TxDb and extracting introns...")
gtf_raw <- rtracklayer::import(stringtie_gtf)

gtf_clean <- gtf_raw[strand(gtf_raw) %in% c("+", "-")]

txdb_st <- txdbmaker::makeTxDbFromGRanges(gtf_clean)
st_introns_grl <- GenomicFeatures::intronsByTranscript(txdb_st, use.names=TRUE)
st_introns_flat <- unlist(st_introns_grl)

# 3. FILTERING
# 3.1. Remove novel events not validated by ds_detection outputs or references
matches_safe <- findOverlaps(st_introns_flat, safe_ranges, type="equal", maxgap=1)
unsupported_indices <- setdiff(seq_along(st_introns_flat), queryHits(matches_safe))
unsupported_tx_names <- unique(names(st_introns_flat[unsupported_indices]))

# 3.2. Remove events without a single novel event (already in reference)
matches_consensus <- findOverlaps(st_introns_flat, consensus_ranges, type="equal", maxgap=1)
supported_by_lsv_tx_names <- unique(names(st_introns_flat[queryHits(matches_consensus)]))

# COMBINE LOGIC
all_st_tx_names <- names(st_introns_grl)
valid_tx_names <- intersect(supported_by_lsv_tx_names, setdiff(all_st_tx_names, unsupported_tx_names))

lgr$info(sprintf("FILTER: %d transcripts passed both filters.", length(valid_tx_names)))

# 4. EXPORT
lgr$info("EXPORT: Saving filtered novel and augmented GTFs...")
gtf_novel_only <- gtf_clean[mcols(gtf_clean)$transcript_id %in% valid_tx_names & grepl("^MSTRG", mcols(gtf_clean)$transcript_id)]
rtracklayer::export(gtf_novel_only, merged_filtered_gtf, format="gtf")

# For the augmented version, we still keep ALL reference transcripts,
# plus our "Highly Validated" novel ones.
ref_pattern <- ifelse(species == "Mus_musculus", "ENSMUST", "ENST")

gtf_with_ref <- c(
    gtf_clean[grepl(ref_pattern, as.character(mcols(gtf_clean)$transcript_id))],
    gtf_novel_only
)

rtracklayer::export(gtf_with_ref, merged_fil_withRef_gtf, format="gtf")
lgr$info("CONSTRUCTING GTF: Complete.")