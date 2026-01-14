# this script collapses all novel transcripts in a gene into a meta-novel transcript
# Output: tx2tx file, tx2gene file that you can use with sleuth

# SNAKEMAKE PARAMS
gtf_merged <- snakemake@input[["gtf_merged"]]
gtf_filtered <- snakemake@input[["gtf_filtered"]]

module_path <- snakemake@params[["module_path"]]
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
# note: do i still need this mapping part?
# import stringtie merged GTF - contains the mappings you need for gene to MSTRG identifier
merged_stringtie_gtf_full <- rtracklayer::import(paste(module_path,gtf_merged,sep=""))
# contains ONLY transcripts with HIGH CONFIDENCE splicing events
filtered_stringtie_gtf <- rtracklayer::import(paste(module_path,gtf_filtered,sep=""))

# load all as a dataframe
merged_stringtie_gtf_full_df <- as.data.frame(merged_stringtie_gtf_full)
filtered_stringtie_gtf_df <- as.data.frame(filtered_stringtie_gtf)

lgr$info("Getting transcript-to-gene mappings in stringtie GTF...")
# get MSTRG:ENSMUSG mappings
mappings <- merged_stringtie_gtf_full_df[,c("gene_id","ref_gene_id")]
mappings <- unique(na.omit(mappings))

lgr$info("Annotating high-confidence filtered transcripts with gene...")
# map the filtered novel transcripts to give them the ENSMUSG column
filtered_stringtie_gtf_df_with_refgeneid <- dplyr::inner_join(filtered_stringtie_gtf_df,mappings,by="gene_id")

# now we have the gene id for the novel transcripts - can combine with t2g
t2g_augment <- filtered_stringtie_gtf_df_with_refgeneid[,c("transcript_id","ref_gene_id")]
t2g_augment <- unique(t2g_augment)
colnames(t2g_augment) <- c("target_id","ens_gene")

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
write.csv(t2g_augment, file=paste(module_path,out_uncollapsed,sep=""), row.names=FALSE)
lgr$info("Uncollapsed t2g dataframe saved at: %s", paste(module_path,out_uncollapsed,sep=""))

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

write.csv(t2g_augment, file=paste(module_path,out_collapsed,sep=""), row.names=FALSE)

lgr$info("Collapsed t2g dataframe saved at: %s", paste(module_path,out_collapsed,sep=""))
lgr$remove_appender("snakemake_file_log")