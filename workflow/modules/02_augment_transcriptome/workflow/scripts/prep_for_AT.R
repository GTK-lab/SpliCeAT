consensus_events <- snakemake@input[["consensus_events_LSVs"]]
stringtie_gtf    <- snakemake@input[["stringtie_gtf"]]

species                <- snakemake@params[["species"]]
merged_filtered_gtf    <- snakemake@output[["merged_filtered_gtf"]]
merged_fil_withRef_gtf <- snakemake@output[["merged_fil_withRef_gtf"]]

log_file <- snakemake@log[[1]]

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

# GTF FILTERING
lgr$info("CONSTRUCTING GTF: Processing merged Stringtie GTF file...")
gtf_raw <- rtracklayer::import(stringtie_gtf)

mcols(gtf_raw)$transcript_id <- gsub("transcript:", "", as.character(mcols(gtf_raw)$transcript_id))
mcols(gtf_raw)$gene_id       <- gsub("gene:", "", as.character(mcols(gtf_raw)$gene_id))

# Keep only features with valid strands for TxDb construction
gtf_clean <- gtf_raw[strand(gtf_raw) %in% c("+", "-")]

lgr$info("CONSTRUCTING GTF: Building Transcript Database from Stringtie GTF...")
txdb             <- txdbmaker::makeTxDbFromGRanges(gtf_clean)
gtf_introns_grl  <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)
gtf_introns_flat <- unlist(gtf_introns_grl)

lgr$info("CONSTRUCTING GTF: Matching Consensus Events to Introns...")
consensus_matrix <- fread(consensus_events) %>% as.data.frame()

consensus_gr <- makeGRangesFromDataFrame(
  consensus_matrix,
  keep.extra.columns = TRUE,
  seqnames.field     = "chr",
  start.field        = "start",
  end.field          = "end"
)

# Find introns that match our validated junctions exactly
hits <- findOverlaps(gtf_introns_flat, consensus_gr, type="equal", maxgap=1)
overlapping_tx_ids <- unique(names(gtf_introns_flat[queryHits(hits)]))

# Extract novel 'MSTRG' transcripts that passed the consensus check
novel_supported_tx <- overlapping_tx_ids[grepl("MSTRG", overlapping_tx_ids)]

lgr$info("CONSTRUCTING GTF: Isolating novel isoforms...")
gtf_novel_only <- gtf_clean[mcols(gtf_clean)$transcript_id %in% novel_supported_tx]
rtracklayer::export(gtf_novel_only, merged_filtered_gtf, format="gtf")

lgr$info(sprintf("CONSTRUCTING GTF: Filtering complete. %d novel transcripts validated.", length(novel_supported_tx)))

# Augmented to reference GTF
ref_pattern <- ifelse(species == "Mus_musculus", "ENSMUST", "ENST")

lgr$info("CONSTRUCTING GTF: Merging novel supported transcripts with reference...")

gtf_with_ref <- c(
  gtf_clean[grepl(ref_pattern, as.character(mcols(gtf_clean)$transcript_id))],
  gtf_novel_only
)

rtracklayer::export(gtf_with_ref, merged_fil_withRef_gtf, format="gtf")
lgr$info(sprintf("CONSTRUCTING GTF: Complete Augmented GTF generated at %s", merged_fil_withRef_gtf))