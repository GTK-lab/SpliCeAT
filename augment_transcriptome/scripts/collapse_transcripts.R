# this script collapses all novel transcripts in a gene into a meta-novel transcript
# Output: tx2tx file, tx2gene file that you can use with sleuth

#### CHANGE THIS ####
config_file_path <- "/mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/config/config.yaml"
organism <- "mouse" # change to "human" if needed
####################

library(lgr)

lgr$info("Loading libraries...")

suppressMessages({
library(yaml)
library(dplyr)
library(rtracklayer)
library(biomaRt)
})

lgr$info("Done.")

# load in config containing file paths and other params
lgr$info("Reading config.yaml...")
config <- read_yaml(config_file_path)
lgr$info("Done.")

lgr$info("Loading in stringtie merged GTF...")
# note: do i still need this mapping part?
# import stringtie merged GTF - contains the mappings you need for gene to MSTRG identifier
merged_stringtie_gtf_full <- rtracklayer::import(paste(config$BASE_PATH,"results/merged_assembly/merged_stringtie_assembly.gtf",sep=""))
# contains ONLY transcripts with HIGH CONFIDENCE splicing events
filtered_stringtie_gtf <- rtracklayer::import(paste(config$BASE_PATH,"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf",sep=""))

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
if (organism == "mouse"){
  ensembl <- useEnsembl(biomart = 'genes', 
                         dataset = 'mmusculus_gene_ensembl',
                         version = config$mouse_ensembl_version)
} else if (organism == "human") {
  ensembl <- useEnsembl(biomart = "genes", 
                        dataset = "hsapiens_gene_ensembl",
                        version = config$hsapiens_ensembl_version)
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
write.csv(t2g_augment, file=paste(config$BASE_PATH,"results/augmented_transcriptome/t2g_augment_uncollapsed.csv",sep=""), row.names=FALSE)
lgr$info("Uncollapsed t2g dataframe saved at: %s", paste(config$BASE_PATH,"results/augmented_transcriptome/t2g_augment_uncollapsed.csv",sep=""))

# function for collaspsing transcripts to create tx to tx group
if (organism == "mouse"){
  collapse_transcripts <- function(row){
  if (grepl("ENSMUST",row[["target_id"]])){
    row[["target_id"]]
  } else {
    s <- strsplit(row[["target_id"]], ".", fixed = TRUE)[[1]]
    paste(s[[1]],".",s[[2]],".","NovelGroup",sep="")
  }
}
  } else if (organism == "human") {
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

write.csv(t2g_augment, file=paste(config$BASE_PATH,"results/augmented_transcriptome/t2g_augment_collapsed.csv",sep=""), row.names=FALSE)
lgr$info("Collapsed t2g dataframe saved at: %s", paste(config$BASE_PATH,"results/augmented_transcriptome/t2g_augment_collapsed.csv",sep=""))
