# Snakemake Parameters
consensus_events <- snakemake@input[["consensus_events_LSVs"]]
stringtie_gtf    <- snakemake@input[["stringtie_gtf"]]
reference_gtf    <- snakemake@input[["reference_gtf"]]

species                <- snakemake@params[["species"]]
novel_gtf    <- snakemake@output[["novel_gtf"]]
augmented_gtf          <- snakemake@output[["augmented_gtf"]]
consensus_GTF_LSVs     <- snakemake@output[["validated_lsv_tsv"]]

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

# 1. TRUTH SET (Ref + ds_detection consensus)
lgr$info("TRUTH SET: Processing Reference GTF and Granges Consensus File...")

official_ref <- rtracklayer::import(reference_gtf)
mcols(official_ref)$transcript_id <- clean_ids(mcols(official_ref)$transcript_id)
mcols(official_ref)$gene_id <- clean_ids(mcols(official_ref)$gene_id)

ref_txdb <- txdbmaker::makeTxDbFromGRanges(official_ref)
ref_introns <- GenomicFeatures::intronsByTranscript(ref_txdb) %>% unlist() %>% unique()
ref_exons <- official_ref[official_ref$type == "exon"]


consensus_matrix <- fread(consensus_events) %>% as.data.frame()

lgr$info("TRUTH SET: Collapsing feature types from Granges Consensus File...")


use_granges <- snakemake@params[["overlap_type"]]

consensus_matrix <- fread(consensus_events) %>% as.data.frame()

if (use_granges) {
  # Event Level Consensus
  lgr$info("TRUTH SET: Collapsing feature types from Granges Consensus File...")
  collapse_feature_type <- function(m, w, l) {
    types <- c(as.character(m), as.character(w), as.character(l))
    types <- gsub("\\[|\\]", "", types)
    if (any(grepl("intron_retention|IR_sign_match", types, ignore.case = TRUE), na.rm = TRUE)) {
      return("intron_retention") }
    if (any(grepl("splice_junction|boundary match", types, ignore.case = TRUE), na.rm = TRUE)) {
      return("splice_junction")}
    return(NA_character_)
  }
  consensus_matrix <- consensus_matrix %>%
    rowwise() %>%
    mutate(event_type = collapse_feature_type(
      majiq_feature_type,
      whippet_feature_type,
      leafcutter_feature_type
    )) %>%
    ungroup()
} else {
  lgr$info("TRUTH SET: Using gene-level consensus feature types directly...")
  consensus_matrix <- consensus_matrix %>%
    rename(event_type = feature_type)
}

consensus_matrix <- consensus_matrix %>%
  select(chr, strand, start, end, gene_id, gene_name, event_type) %>%
  filter(!is.na(event_type))


lgr$info("TRUTH SET: Extracting ranges from Granges Consensus File...")
consensus_ranges <- GRanges(
    seqnames = consensus_matrix$chr,
    ranges = IRanges(start = consensus_matrix$start, end = consensus_matrix$end),
    strand = consensus_matrix$strand,
	event_type = consensus_matrix$event_type
)

sj_consensus <- consensus_ranges[consensus_ranges$event_type == "splice_junction"]
lgr$info("TRUTH SET: Splice Junctions Extracted.")
ir_consensus <- consensus_ranges[consensus_ranges$event_type == "intron_retention"]
lgr$info("TRUTH SET: Intron Retention events Extracted.")

exon_consensus <- consensus_ranges[consensus_ranges$event_type == "exon_node"]
if(!use_granges) {
	lgr$info("TRUTH SET: Exon Nodes Extracted.")
}

# 2. TEST SET (StringTie)
lgr$info("STRINGTIE: Building TxDb and extracting genomic features...")
gtf_raw <- rtracklayer::import(stringtie_gtf)
gtf_clean <- gtf_raw[strand(gtf_raw) %in% c("+", "-")]
mcols(gtf_clean)$transcript_id <- clean_ids(mcols(gtf_clean)$transcript_id)
mcols(gtf_clean)$gene_id <- clean_ids(mcols(gtf_clean)$gene_id)

txdb_st <- txdbmaker::makeTxDbFromGRanges(gtf_clean)
st_introns_grl <- GenomicFeatures::intronsByTranscript(txdb_st, use.names=TRUE)
st_introns_flat <- unlist(st_introns_grl)

#extracting exons
st_exons_grl <- GenomicFeatures::exonsBy(txdb_st, by="tx", use.names=TRUE)
st_exons_flat <- unlist(st_exons_grl)


all_st_tx_names <- names(st_introns_grl)
novel_st_tx_names <- all_st_tx_names[grepl("^MSTRG", all_st_tx_names)]


# 3. FILTERING
lgr$info("FILTER: Validating Splice Junctions and Intron Retentions...")

# 3.1. Remove transcripts with stringtie-only novel events
# 3.1.1 Event is splice junction
safe_SJ_ranges <- unique(c(ref_introns, sj_consensus)) #combined truth set
matches_safe <- findOverlaps(st_introns_flat, safe_SJ_ranges, type="equal", maxgap=1)

unsupported_sj_indices <- setdiff(seq_along(st_introns_flat), queryHits(matches_safe)) # stringtie SJ minus truth set SJ
unsupported_sj_tx <- unique(names(st_introns_flat[unsupported_sj_indices]))
lgr$info("FILTER: Splice Junction Validation Complete. Removed %d transcripts with unsupported splice junctions.", length(unsupported_sj_tx))

# 3.1.2 Event is intron retention
# case 1: intron retention in stringtie is otherwise an intron in reference
matches_ref_skipping <- findOverlaps(st_exons_grl, ref_introns, type="within")
tx_skipping_ref_introns <- unique(names(st_exons_grl[queryHits(matches_ref_skipping)]))

# case 2: intron retention in stringtie matched DS output
matches_ir_validation <- findOverlaps(st_exons_grl, ir_consensus, type="within")
tx_with_ir_support <- unique(names(st_exons_grl[queryHits(matches_ir_validation)]))

unsupported_ir_tx <- setdiff(tx_skipping_ref_introns, tx_with_ir_support) #case 1 minus case 2
lgr$info("FILTER: Intron Retention Validation Complete. Removed %d transcripts with unsupported intron retention events.", length(unsupported_ir_tx))

# 3.1.3 Event is exon node (for gene_overlap input)
tx_with_exon_support <- character(0)

if (!use_granges) {
    lgr$info("FILTER: Validating Exon Nodes")

	# identify transcripts with at least one supported exon node
    matches_exon_lsv <- findOverlaps(st_exons_grl, exon_consensus, type="equal")
    tx_with_exon_support <- unique(names(st_exons_grl[queryHits(matches_exon_lsv)]))
}

# COMBINED UNSUPPORTED TRANSCRIPTS
unsupported_tx_names <- unique(c(unsupported_sj_tx, unsupported_ir_tx))
lgr$info("FILTER: %d transcripts contains unvalidated events. They have been removed", length(unsupported_tx_names))

# 3.2. Remove events without a single novel event (already in reference)
lgr$info("FILTER: Checking for novel events...")
# 3.2.1 Event is splice junction
matches_sj_lsv <- findOverlaps(st_introns_flat, sj_consensus, type="equal", maxgap=1)
tx_with_sj_support <- unique(names(st_introns_flat[queryHits(matches_sj_lsv)]))

# 3.2.2 Event is intron retention
matches_ir_lsv <- findOverlaps(st_exons_grl, ir_consensus, type="within") #intron spans within long exon
tx_with_ir_support <- unique(names(st_exons_grl[queryHits(matches_ir_lsv)]))

supported_by_truth <- unique(c(tx_with_sj_support, tx_with_ir_support, tx_with_exon_support))

# 3.3 Combined filters to identify valid transcripts
valid_tx_names <- intersect(supported_by_truth, setdiff(novel_st_tx_names, unsupported_tx_names))
lgr$info(sprintf("FILTER: Overall Filtering Complete. %d transcripts contains at least one validated novel event.", length(valid_tx_names)))


# 4. EXPORT
# 4.1. EXPORT CONSENSUS MATRIX
lgr$info("EXPORT: Generating validated consensus matrix...")

#Subset features to the final "winning" novel transcripts
st_final_introns <- st_introns_flat[names(st_introns_flat) %in% valid_tx_names]
st_final_exons <- st_exons_grl[names(st_exons_grl) %in% valid_tx_names]

# Match consensus events back to these novel transcripts
valid_sj_idx <- unique(queryHits(findOverlaps(sj_consensus, st_final_introns, type="equal", maxgap=1)))
valid_ir_idx <- unique(queryHits(findOverlaps(ir_consensus, st_final_exons, type="within")))

valid_exon_idx <- if(!use_granges) {
    unique(queryHits(findOverlaps(exon_consensus, st_final_exons, type="equal")))
} else { numeric(0) }

final_consensus_matrix <- rbind(
    consensus_matrix %>% filter(event_type == "splice_junction") %>% slice(valid_sj_idx),
    consensus_matrix %>% filter(event_type == "intron_retention") %>% slice(valid_ir_idx),
    consensus_matrix %>% filter(event_type == "exon_node") %>% slice(valid_exon_idx)
)

fwrite(final_consensus_matrix,consensus_GTF_LSVs, sep="\t")
lgr$info(sprintf("EXPORT: Consensus Matrix Complete. File Saved at %s",consensus_GTF_LSVs))

# 4.2. EXPORT GTF
# 4.2.1. Novel transcripts only (filtered StringTie)
lgr$info("EXPORT: Generating validated consensus matrix (supported by StringTie)...")
gtf_novel_only <- gtf_clean[mcols(gtf_clean)$transcript_id %in% valid_tx_names]
rtracklayer::export(gtf_novel_only, novel_gtf, format="gtf")
lgr$info(sprintf("EXPORT: Novel-only GTF Construction Complete. File Saved at %s",novel_gtf))

# 4.2.2. Augmented GTF (Ref + Validated Novel)
# ref_pattern <- ifelse(species == "Mus_musculus", "ENSMUST", "ENST")
gtf_with_ref <- c(
    official_ref,
	gtf_novel_only
)
rtracklayer::export(gtf_with_ref, augmented_gtf, format="gtf")
lgr$info(sprintf("EXPORT: Augmented GTF Construction Complete. File Saved at %s", augmented_gtf))