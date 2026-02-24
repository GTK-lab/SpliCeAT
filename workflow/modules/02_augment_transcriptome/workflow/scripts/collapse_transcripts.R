# SNAKEMAKE PARAMS
gtf_merged   <- snakemake@input[["gtf_merged"]]
gtf_filtered <- snakemake@input[["gtf_filtered"]]
organism     <- snakemake@params[["organism"]]
ensembl_ver  <- snakemake@params[["ensembl"]] # Renamed to avoid conflict with biomaRt object
out_uncollapsed <- snakemake@output[["uncollapsed"]]
out_collapsed   <- snakemake@output[["collapsed"]]

log_file <- snakemake@log[[1]]

# LIBRARIES
suppressMessages({
  library(dplyr)
  library(rtracklayer)
  library(biomaRt)
  library(data.table)
  library(lgr)
})

library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")

lgr$info("Loading and cleaning GTFs...")

# Helper to strip prefixes and version decimals
clean_ids <- function(x) {
  x <- gsub("transcript:|gene:", "", as.character(x))
  x <- ifelse(grepl("^MSTRG", x),
              # For MSTRG: Keep the FULL ID (MSTRG.123.1) so it stays unique for gffread
              x,
              # For Ensembl: Strip the version (ENSMUST000.14 -> ENSMUST000)
              sub("\\..*$", "", x)
  )
  return(x)
}

stringtie_gtf_df <- as.data.frame(rtracklayer::import(gtf_merged)) %>%
  mutate(across(c(transcript_id, gene_id, ref_gene_id), clean_ids))

filtered_gtf_df <- as.data.frame(rtracklayer::import(gtf_filtered)) %>%
  mutate(across(c(transcript_id, gene_id), clean_ids))

lgr$info("Mapping MSTRG to Ensembl Gene IDs...")
# Create a map between StringTie gene IDs and Reference Gene IDs
gene_map <- stringtie_gtf_df %>%
  filter(!is.na(ref_gene_id) & ref_gene_id != "") %>%
  dplyr::select(gene_id, ref_gene_id) %>%
  distinct()

lgr$info("Fetching Reference annotations from BioMart...")
mart <- useEnsembl(biomart = "genes",
                   dataset = ifelse(organism == "Mus_musculus", "mmusculus_gene_ensembl", "hsapiens_gene_ensembl"),
                   version = ensembl_ver)

# Get reference t2g (without versions)
t2g_ref <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart) %>%
  dplyr::rename(target_id = ensembl_transcript_id,
                ens_gene = ensembl_gene_id,
                ext_gene = external_gene_name)

lgr$info("Building Augmented T2G...")
# Prepare novel transcripts
t2g_novel <- filtered_gtf_df %>%
  filter(type == "transcript" & grepl("MSTRG", transcript_id)) %>%
  dplyr::select(target_id = transcript_id, gene_id) %>%
  distinct() %>%
  left_join(gene_map, by = "gene_id") %>%
  # If no ref_gene_id found, keep the MSTRG gene ID
  mutate(ens_gene = ifelse(is.na(ref_gene_id), gene_id, ref_gene_id)) %>%
  # Join with reference gene names (Left join ensures orphans aren't dropped)
  left_join(distinct(t2g_ref, ens_gene, ext_gene), by = "ens_gene") %>%
  mutate(ext_gene = ifelse(is.na(ext_gene), ens_gene, ext_gene)) %>%
  dplyr::select(target_id, ens_gene, ext_gene)

# Combine Ref and Novel
t2g_final <- bind_rows(t2g_ref, t2g_novel) %>% distinct(target_id, .keep_all = TRUE)

# Fill empty gene names with the Gene ID so your plots aren't blank
t2g_final <- t2g_final %>%
  mutate(ext_gene = ifelse(is.na(ext_gene) | ext_gene == "", ens_gene, ext_gene))

write.csv(t2g_final, file = out_uncollapsed, row.names = FALSE)
lgr$info("Uncollapsed T2G saved.")

lgr$info("Creating Collapsed Groups...")
# Vectorized collapse logic
t2g_final <- t2g_final %>%
  mutate(collapsed_target_id = ifelse(grepl("MSTRG", target_id),
                                      paste0(ens_gene, ".NovelGroup"),
                                      target_id))

write.csv(t2g_final, file = out_collapsed, row.names = FALSE)
lgr$info("Collapsed T2G saved.")