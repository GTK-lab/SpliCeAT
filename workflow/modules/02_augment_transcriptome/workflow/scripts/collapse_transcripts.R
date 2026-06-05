# SNAKEMAKE PARAMS
gtf_stringtie   <- snakemake@input[["gtf_stringtie"]]
gtf_novel_only <- snakemake@input[["gtf_novel_only"]]
organism     <- snakemake@params[["organism"]]
ensembl_ver  <- snakemake@params[["ensembl"]] # Renamed to avoid conflict with biomaRt object
out_uncollapsed <- snakemake@output[["uncollapsed"]]
out_collapsed   <- snakemake@output[["collapsed"]]

log_file <- snakemake@log[[1]]
library(lgr)
log_con <- file(log_file, open = "a")
sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")

# LIBRARIES
suppressMessages({
  library(dplyr)
  library(rtracklayer)
  library(biomaRt)
  library(data.table)
  library(lgr)
})

lgr$info("LOADING: Retrieving and cleaning GTFs...")
# Cleaning Metadata columns
clean_ids <- function(x) {
	x <- gsub("transcript:|gene:", "", as.character(x))
	x <- ifelse(grepl("^MSTRG", x), x, sub("\\..*$", "", x))
	return(x)}

# Load GTFs as dataframe + cleaned IDs
stringtie_gtf_df <- as.data.frame(rtracklayer::import(gtf_stringtie)) %>%
  mutate(across(c(transcript_id, gene_id, ref_gene_id), clean_ids))

filtered_gtf_df <- as.data.frame(rtracklayer::import(gtf_novel_only)) %>%
  mutate(across(c(transcript_id, gene_id), clean_ids))
lgr$info("LOADING: Complete.")

lgr$info("Fetching Reference annotations from BioMart...")
# Connect to Ensembl BioMart
mart <- useEnsembl(biomart = "genes",
                   dataset = ifelse(organism == "Mus_musculus", "mmusculus_gene_ensembl", "hsapiens_gene_ensembl"),
                   version = ensembl_ver)

# Reference Transcript2Gene Map file
t2g_ref <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart) %>%
  dplyr::rename(target_id = ensembl_transcript_id,
                ens_gene = ensembl_gene_id,
                ext_gene = external_gene_name)

lgr$info("Building Augmented T2G...")

# using stringtie gtf to retrieve mapping between MSTRG id and Ensembl ID
lgr$info("Mapping MSTRG IDs to Ensembl Gene IDs...")
gene_map <- stringtie_gtf_df %>%
  filter(!is.na(ref_gene_id) & ref_gene_id != "") %>%
  dplyr::select(gene_id, ref_gene_id) %>%
  distinct()

# using final filtered gtf to retrieve mapping between transcript id and gene id
t2g_novel <- filtered_gtf_df %>%
  filter(type == "transcript" & grepl("MSTRG", transcript_id)) %>%
  dplyr::select(target_id = transcript_id, gene_id) %>%
  distinct() %>%

  # cross-reference with gene map for ref_gene_ids, else keep the MSTRG gene ID
  left_join(gene_map, by = "gene_id") %>%
  mutate(ens_gene = ifelse(is.na(ref_gene_id), gene_id, ref_gene_id)) %>%

  # cross reference with biomart for gene name
  left_join(distinct(t2g_ref, ens_gene, ext_gene), by = "ens_gene") %>%
  mutate(ext_gene = ifelse(is.na(ext_gene), ens_gene, ext_gene)) %>%
  dplyr::select(target_id, ens_gene, ext_gene)

# Final: Ref + Novel
t2g_final <- bind_rows(t2g_ref, t2g_novel) %>%
	# remove duplicates
	distinct(target_id, .keep_all = TRUE) %>%
	# if Gene_Name is NA, add Gene_Id
	mutate(ext_gene = ifelse(is.na(ext_gene) | ext_gene == "", ens_gene, ext_gene))

write.csv(t2g_final, file = out_uncollapsed, row.names = FALSE)
lgr$info(sprintf("Uncollapsed T2G saved at %s",(out_uncollapsed)))

lgr$info("Creating Collapsed Groups...")
# Vectorized collapse logic

t2g_final <- t2g_final %>%
  mutate(collapsed_target_id = ifelse(grepl("MSTRG", target_id),
                                      paste0(ens_gene, ".NovelGroup"),
                                      target_id))

write.csv(t2g_final, file = out_collapsed, row.names = FALSE)

lgr$info("Collapsed T2G saved at %s",(out_collapsed))