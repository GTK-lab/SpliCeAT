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

clean_ids <- function(x) {
  x <- gsub("transcript:|gene:", "", as.character(x))
  x <- ifelse(grepl("^MSTRG", x), x, sub("\\..*$", "", x))
  return(x)
}

mcols(gtf_raw)$transcript_id <- clean_ids(mcols(gtf_raw)$transcript_id)
mcols(gtf_raw)$gene_id       <- clean_ids(mcols(gtf_raw)$gene_id)

if ("ref_gene_id" %in% colnames(mcols(gtf_raw))) {
    mcols(gtf_raw)$ref_gene_id <- clean_ids(mcols(gtf_raw)$ref_gene_id)
}

# Keep only features with valid strands for TxDb construction
gtf_clean <- gtf_raw[strand(gtf_raw) %in% c("+", "-")]

lgr$info("CONSTRUCTING GTF: Building Transcript Database from Stringtie GTF...")
txdb             <- txdbmaker::makeTxDbFromGRanges(gtf_clean)
gtf_introns_grl  <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)
gtf_introns_flat <- unlist(gtf_introns_grl)

gtf_exons <- gtf_clean[gtf_clean$type == "exon"]

lgr$info("CONSTRUCTING GTF: Matching Consensus Events...")
consensus_matrix <- fread(consensus_events) %>% as.data.frame()

# Splice Junctions (start/end match)
junctions_gr <- makeGRangesFromDataFrame(
    consensus_matrix[grepl("junction", consensus_matrix$majiq_type) |
                     grepl("junction", consensus_matrix$whippet_type) |
                     grepl("junction", consensus_matrix$leaf_type), ],
    keep.extra.columns = TRUE, seqnames.field = "chr"
)

hits_introns <- findOverlaps(gtf_introns_flat, junctions_gr, type="equal", maxgap=1)
tx_by_junctions <- unique(names(gtf_introns_flat[queryHits(hits_introns)]))
novel_tx_junc <- tx_by_junctions[grepl("MSTRG", tx_by_junctions)]
lgr$info("CONSTRUCTING GTF: Junction Validation Complete. Filtered %d transcripts", length(novel_tx_junc))


# Nodes/IR (Overlap match)
nodes_gr <- makeGRangesFromDataFrame(
    consensus_matrix[grepl("exon_node|intron_retention", consensus_matrix$majiq_type) |
                     grepl("exon_node|intron_retention", consensus_matrix$whippet_type), ],
    keep.extra.columns = TRUE, seqnames.field = "chr"
)

hits_exons <- findOverlaps(gtf_exons, nodes_gr, type="any")
tx_by_nodes <- unique(mcols(gtf_exons[queryHits(hits_exons)])$transcript_id)
novel_tx_node <- tx_by_nodes[grepl("MSTRG", tx_by_nodes)]
lgr$info("CONSTRUCTING GTF: Node Validation Complete. Filtered %d transcripts", length(novel_tx_node))

# Find introns that match our validated junctions exactly
overlapping_tx_ids <- unique(c(tx_by_junctions, tx_by_nodes))
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