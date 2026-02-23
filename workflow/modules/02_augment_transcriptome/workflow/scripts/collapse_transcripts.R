# this script collapses all novel transcripts in a gene into a meta-novel transcript
# Output: tx2tx file, tx2gene file that you can use with sleuth

# SNAKEMAKE PARAMS
gtf_merged <- snakemake@input[["gtf_merged"]]
gtf_filtered <- snakemake@input[["gtf_filtered"]]

organism <- snakemake@params[["organism"]]
ensembl <- snakemake@params[["ensembl"]]

out_uncollapsed <- snakemake@output[["uncollapsed"]]
out_collapsed <- snakemake@output[["collapsed"]]
n_threads <- snakemake@threads

# LOGGER SETUP
library(lgr)
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "a")

sink(log_con, append = FALSE)
sink(log_con, append = FALSE, type = "message")

lgr$set_appenders(list())
lgr$set_propagate(FALSE)

lgr$add_appender(AppenderFile$new(log_file, threshold = "all"), name = "snakemake_file_log")

options(lgr.log_messages = TRUE)
options(lgr.log_warnings = TRUE)
options(warn = 1)

# LIBRARIES
lgr$info("Loading libraries...")

suppressMessages({
library(yaml)
library(dplyr)
library(rtracklayer)
library(biomaRt)
})

lgr$info("Done.")

lgr$info("Loading in stringtie merged GTF...")
# load all as a dataframe
stringtie_gtf_df <- as.data.frame(rtracklayer::import(paste(gtf_merged))) #stringtie final output for MSTRG identifier
filtered_gtf_df <- as.data.frame(rtracklayer::import(paste(gtf_filtered))) #augmented transcriptome -> only novel high-confidence transcripts

# cleaning for cross reference
clean_ids <- function(x) {
  x <- gsub("transcript:", "", as.character(x))
  x <- gsub("gene:", "", as.character(x))
  return(x)
}

stringtie_gtf_df$transcript_id <- clean_ids(stringtie_gtf_df$transcript_id)
stringtie_gtf_df$gene_id       <- clean_ids(stringtie_gtf_df$gene_id)
stringtie_gtf_df$ref_gene_id   <- clean_ids(stringtie_gtf_df$ref_gene_id)

filtered_gtf_df$transcript_id    <- clean_ids(filtered_gtf_df$transcript_id)
filtered_gtf_df$gene_id          <- clean_ids(filtered_gtf_df$gene_id)


lgr$info("Getting transcript-to-gene mappings in stringtie GTF...")
# get MSTRG:ENSMUSG mappings
mappings <- stringtie_gtf_df %>%
  dplyr::select(gene_id, ref_gene_id) %>%
  distinct()


lgr$info("Annotating high-confidence filtered transcripts with gene...")
t2g_augment <- filtered_gtf_df %>%
  filter(type == "transcript") %>%
  dplyr::select(target_id = transcript_id, gene_id) %>%
  distinct() %>%
  left_join(mappings, by = "gene_id") %>%
  # If there is no ref_gene_id, use the MSTRG gene_id as a placeholder
  mutate(ens_gene = ifelse(is.na(ref_gene_id) | ref_gene_id == "", gene_id, ref_gene_id)) %>%
  dplyr::select(target_id, ens_gene)

lgr$info(sprintf("Processed %d novel transcripts for augmentation.", nrow(t2g_novel)))

lgr$info("Getting ensembl annotations...")

# load in normal t2g
if (organism == "Mus_musculus"){
  ensembl <- useEnsembl(biomart = 'genes',
                         dataset = 'mmusculus_gene_ensembl',
                         version = ensembl)
} else if (organism == "Homo_sapiens"){
  ensembl <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = ensembl)
}

t2g <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version",
                            "external_gene_name"),mart = ensembl)

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id_version,
                     ens_gene = ensembl_gene_id_version, ext_gene = external_gene_name)

# mapping of ens_gene to ext_gene (ENSMUSG:gene name)
ens_gene_ext_gene <- unique(t2g[,c(2,3)])

# augment my high-confidence transcripts with gene name
t2g_augment <- dplyr::inner_join(t2g_augment,ens_gene_ext_gene,by="ens_gene")

lgr$info("Combining reference transcripts and high-confidence novel transcripts into a t2g dataframe...")
# now we rbind it with normal t2g to create an AUGMENTED T2G
t2g_augment <- rbind(t2g, t2g_augment)
head(t2g_augment)
tail(t2g_augment)
write.csv(t2g_augment, file=paste(out_uncollapsed), row.names=FALSE)
lgr$info("Uncollapsed t2g dataframe saved at: %s", paste(out_uncollapsed,sep=""))

# function for collaspsing transcripts to create tx to tx group
if (organism == "Mus_musculus") {
  collapse_transcripts <- function(row){
  if (grepl("ENSMUST",row[["target_id"]])){
    row[["target_id"]]
  } else {
    s <- strsplit(row[["target_id"]], ".", fixed = TRUE)[[1]]
    paste(s[[1]],".",s[[2]],".","NovelGroup",sep="")
  }
}
  } else if (organism == "Homo_sapiens") {
  collapse_transcripts <- function(row){
  if (grepl("ENST",row[["target_id"]])){
    row[["target_id"]]
  } else {
    s <- strsplit(row[["target_id"]], ".", fixed = TRUE)[[1]]
    paste(s[[1]],".",s[[2]],".","NovelGroup",sep="")
  }
}
  }

t2g_augment$collapsed_target_id <- apply(t2g_augment,1,collapse_transcripts)
# remove duplicated target_id entries becos there can be overlapping genes
t2g_augment <- t2g_augment[!duplicated(t2g_augment[c('target_id')]),]

write.csv(t2g_augment, file=paste(out_collapsed), row.names=FALSE)

lgr$info("Collapsed t2g dataframe saved at: %s", paste(out_collapsed))
lgr$remove_appender("snakemake_file_log")