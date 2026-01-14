### prep_for_AT.R
# Prepares splicing masterlists
# Filters merged Stringtie GTF with high confidence splice events

# SNAKEMAKE PARAMS
ds_files <- snakemake@input[["ds_files"]]
majiq_het_tsv <- ds_files[1]
whippet_delta_psi <- ds_files[2]
leafcutter_cluster_file <- ds_files[3]
leafcutter_effect_size_file <- ds_files[4]

gtf_merged <- snakemake@input[["gtf_merged"]]

merged_filtered_gtf <- snakemake@output[["merged_filtered_gtf"]]
merged_fil_withRef_gtf <- snakemake@output[["merged_fil_withRef_gtf"]]

module_path <- snakemake@params[["module_path"]]
species    <- snakemake@params[["organism"]]
version <- snakemake@params[["ensembl"]]
n_threads   <- snakemake@threads

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
library(biomaRt)
library(stringr)
library(data.table)
library(VennDetail)
suppressMessages(library(dplyr))
library(GenomicRanges)
library(rtracklayer)
library(tidyr)
})

lgr$info("Done.")

lgr$info("MASTERLISTS PREP")

### ENSEMBL ANNOTATIONS
lgr$info("Getting Ensembl annotations...")

target_dataset <- switch(
  species,
  Mus_musculus = "mmusculus_gene_ensembl",
  Homo_sapiens = "hsapiens_gene_ensembl",
  stop(paste0("Unsupported organism: %s", species))
)

ensembl <- NULL
try({
  ensembl <- useEnsembl(
    biomart = "genes",
    dataset = target_dataset,
    version = version
  )
}, silent = TRUE)


if (is.null(ensembl)) {
  host <- "https://www.ensembl.org"

  marts <- listMarts(host = host)
  mart_name <- marts$biomart[grep("ENSEMBL_MART_ENSEMBL", marts$biomart)][1]
  datasets <- listDatasets(useMart(mart_name, host = host))

  if (!(target_dataset %in% datasets$dataset)) {
    stop(paste0("Dataset '", target_dataset, "' not found in live Ensembl."))
  }

  ensembl <- useMart(
    biomart = mart_name,
    dataset = target_dataset,
    host = host
  )
}

# Inform which Ensembl was used
if (exists("host")) {
	lgr$info("Connected to live Ensembl at %s", paste(host))
} else {
	lgr$info("Connected to versioned Ensembl: %s release ", paste(version))
}

# Retrieve annotations from the mart
annotations <- getBM(
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name',
    'description',
    'chromosome_name',
    'start_position',
    'end_position',
    'strand'
  ),
  mart = ensembl
)
lgr$info("Done.")

### MAJIQ MASTERLIST
lgr$info("Processing Majiq LSV file...")

# split lsv columns based on colon and select necessary columns
majiq_ids <- read.table(majiq_het_tsv, sep = "\t",header=TRUE, stringsAsFactors = FALSE)$lsv_id
majiq_lsv_final <- data.frame(lsv_id = majiq_ids, stringsAsFactors = FALSE) |>
    tidyr::separate(
        lsv_id,
        into = c("prefix", "Ensembl_ID", "type", "Coord"),
        sep = ":",
        remove = TRUE
    )

majiq_lsv_final <- majiq_lsv_final |>
    mutate(
        gene_ids = str_replace(Ensembl_ID, pattern = "\\.[0-9]+$", replacement = "")
    ) |>
    dplyr::select(gene_ids, Ensembl_ID, type, Coord, prefix)

# write to majiq_masterlist
dir.create(paste(module_path,"masterlists",sep=""))
write.csv(majiq_lsv_final,paste(module_path,"masterlists/majiq_masterlist.csv",sep=""),row.names = FALSE)
lgr$info("Majiq masterlist created at: %s",paste(module_path,"masterlists/majiq_masterlist.csv",sep=""))


### LEAFCUTTER MASTERLIST
lgr$info("Processing Leafcutter input files...")
leafcutter_cluster_sig <- read.table(leafcutter_cluster_file,sep="\t",header=TRUE)
leafcutter_cluster_sig <- leafcutter_cluster_sig[,c(1,6,7)]

leafcutter_effect_sizes <- read.table(leafcutter_effect_size_file,sep="\t",header=TRUE)
leafcutter_effect_sizes[c("chr","start","stop","cluster")] <- str_split_fixed(leafcutter_effect_sizes$intron, ':', 4)

leafcutter_effect_sizes$start <- as.numeric(leafcutter_effect_sizes$start)
leafcutter_effect_sizes$stop <- as.numeric(leafcutter_effect_sizes$stop)

leafcutter_effect_sizes$cluster_joined <- paste(leafcutter_effect_sizes$chr,":",leafcutter_effect_sizes$cluster,sep="")
leafcutter_effect_sizes <- leafcutter_effect_sizes[,c(1,5:10)]

leafcutter_final <- left_join(leafcutter_effect_sizes,
                              leafcutter_cluster_sig,
                              by = c("cluster_joined"="cluster"))

leafcutter_final <- left_join(leafcutter_final, annotations, by = c("genes"="ensembl_gene_id"))

leafcutter_final <- leafcutter_final |> dplyr::rename(intron_start = start, intron_end = stop, Ensembl_ID = genes, P_Adjust = p.adjust, Genes = external_gene_name ) |>
    dplyr::select(
    Ensembl_ID,
    chr,
    Genes,
    intron_start,
    intron_end,
    deltapsi,
    P_Adjust,
    intron
  )

leafcutter_final_filtered <- filter(leafcutter_final, P_Adjust <0.05 & abs(deltapsi) >= 0.2)

write.csv(leafcutter_final_filtered,paste(module_path,"masterlists/leafcutter_masterlist.csv",sep=""),row.names = FALSE)
lgr$info("Leafcutter masterlist created at: %s",paste(module_path,"masterlists/leafcutter_masterlist.csv",sep=""))

### WHIPPET MASTERLIST
lgr$info("Processing Whippet input file...")

# read whippet diff file and select necessary columns
whippet_diff <- read.table(whippet_delta_psi, sep = '\t', header = TRUE, row.names=NULL)
colnames(whippet_diff) <- c("Gene", "Node", "Coord", "Strand", "Type", "Psi_A", "Psi_B", "DeltaPsi", "Probability", "Complexity", "Entropy")
whippet_diff <- whippet_diff[, c("Gene","Coord","Type","DeltaPsi","Probability")]

# clean Gene IDs to match Ensembl ID format
whippet_diff <- whippet_diff |>
  mutate(
    gene_ids = sub("^gene:", "", Gene),
    gene_ids = sub("\\.[0-9]+$", "", gene_ids)
  )

whippet_diff_dt <- as.data.table(whippet_diff)
whippet_diff_highest_prob <- whippet_diff_dt[
  whippet_diff_dt[, .I[Probability == max(Probability)], by=gene_ids]$V1
]

# obtain chromosome ranges
whippet_diff_highest_prob <- whippet_diff_highest_prob |>
  tidyr::separate(Coord, into = c("CE_chr","CE_range"), sep = ":",remove = FALSE) |>
  tidyr::separate(CE_range, into = c("CE_start","CE_stop"), sep = "-")

whippet_diff_final <- whippet_diff_highest_prob %>%
  left_join(
    annotations,
    by = c("gene_ids" = "ensembl_gene_id")
  )

whippet_diff_final <- whippet_diff_final |>
  dplyr::rename(
    Ensembl_ID = gene_ids,
    Ensembl_ID_version = Gene,
    CE_coord = Coord,
    Type = Type,
    DeltaPSI = DeltaPsi,
    Probability = Probability,
    CE_chr = CE_chr,
    CE_start = CE_start,
    CE_stop = CE_stop,
    Gene_name = external_gene_name,
    Description = description,
    Gene_chr = chromosome_name,
    Gene_start = start_position,
    Gene_stop = end_position,
    Gene_strand = strand
  )

# replace -1 and 1 with - and +
whippet_diff_final$Gene_strand[whippet_diff_final$Gene_strand=="-1"]<-"-"
whippet_diff_final$Gene_strand[whippet_diff_final$Gene_strand=="1"]<-"+"

write.csv(whippet_diff_final, paste(module_path,"masterlists/whippet_masterlist.csv",sep=""),row.names = FALSE)

lgr$info("Whippet masterlist created at: %s",paste(module_path,"masterlists/whippet_masterlist.csv",sep=""))

### GET INTERSECTIONS
# get intersections of majiq, whippet, leafcutter
lgr$info("Getting intersections of Majiq, Whippet, Leafcutter...")

# get all Gene Names that pass filtering from each masterlist
filter_genes <- function(whippet_csv, majiq_csv, leafcutter_csv){
  whippet_df <- read.csv(whippet_csv,header=TRUE,check.names=FALSE)
  majiq_df <- read.csv(majiq_csv,header=TRUE,check.names=FALSE)
  leafcutter_df <- read.csv(leafcutter_csv,header=TRUE,check.names=FALSE)

  whippet_df_filtered <- filter(whippet_df, abs(DeltaPSI) >= 0.2 & Probability >= 0.95)
  whippet_df_filtered <- whippet_df_filtered[!(is.na(whippet_df_filtered$Gene_name) | whippet_df_filtered$Gene_name==""), ]
  whippet_genes <- sort(unique(na.omit(whippet_df_filtered[["Gene_name"]])))

  majiq_df <- majiq_df[!(is.na(majiq_df$external_gene_name) | majiq_df$external_gene_name==""), ]
  majiq_genes <- sort(unique(na.omit(majiq_df[["external_gene_name"]])))

  leafcutter_df <- leafcutter_df[!(is.na(leafcutter_df$Genes) | leafcutter_df$Genes==""), ]
  leafcutter_genes <- sort(unique(na.omit(leafcutter_df[["Genes"]])))

  result <- list("majiq_genes" = majiq_genes, "whippet_genes" = whippet_genes, "leafcutter_genes" = leafcutter_genes)
}

filtered_genes <- filter_genes(paste(module_path,"masterlists/whippet_masterlist.csv",sep=""),
                               paste(module_path,"masterlists/majiq_masterlist.csv",sep=""),
                               paste(module_path,"masterlists/leafcutter_masterlist.csv",sep=""))

# create venn object
ven <- venndetail(list(Whippet = filtered_genes$whippet_genes,
                       MAJIQ = filtered_genes$majiq_genes,
                       Leafcutter = filtered_genes$leafcutter_genes))

# Getting elements in subset
shared <- as.data.frame(getSet(ven, subset = c("Shared"))$Detail)
colnames(shared) <- "Gene_name"

MAJIQ_Leafcutter <- as.data.frame(getSet(ven, subset = c("MAJIQ_Leafcutter"))$Detail)
colnames(MAJIQ_Leafcutter) <- "Gene_name"
Whippet_Leafcutter <- as.data.frame(getSet(ven, subset = c("Whippet_Leafcutter"))$Detail)
colnames(Whippet_Leafcutter) <- "Gene_name"
Whippet_MAJIQ <- as.data.frame(getSet(ven, subset = c("Whippet_MAJIQ"))$Detail)
colnames(Whippet_MAJIQ) <- "Gene_name"

MAJIQ <- as.data.frame(getSet(ven, subset = c("MAJIQ"))$Detail)
colnames(MAJIQ) <- "Gene_name"
Leafcutter <- as.data.frame(getSet(ven, subset = c("Leafcutter"))$Detail)
colnames(Leafcutter) <- "Gene_name"
Whippet <- as.data.frame(getSet(ven, subset = c("Whippet"))$Detail)
colnames(Whippet) <- "Gene_name"

# only want union of intersect
union <- rbind(shared, MAJIQ_Leafcutter,Whippet_Leafcutter,Whippet_MAJIQ)
union <- as.data.frame(union[order(union$Gene_name),])
colnames(union) <- "Gene_name"

# merging with whippet masterlist to get event details
master_whippet_df <- read.csv(paste(module_path,"masterlists/whippet_masterlist.csv",sep=""),header=TRUE,check.names=FALSE)
master_whippet_df$Description <- sub(" \\[.*", "", master_whippet_df$Description)

df_to_display <- merge(union, master_whippet_df, by = 'Gene_name')
df_to_display <- df_to_display[,c("Gene_name", "Description","CE_coord", "CE_chr","CE_start","CE_stop","Type","DeltaPSI","Probability", "Gene_strand")]
df_to_display <- df_to_display |>
  dplyr::rename(
    Coordinates = CE_coord,
    Strand = Gene_strand
  )
write.csv(df_to_display, paste(module_path,"masterlists/union_of_intersects_events.csv",sep=""),row.names = FALSE)
lgr$info("Intersections file created at: %s",paste(module_path,"masterlists/union_of_intersects_events.csv",sep=""))

### FILTERING FOR SPLICE EVENTS
lgr$info("FILTERING FOR SPLICE EVENTS")
lgr$info("Processing merged Stringtie GTF file...")

gtf_full <- rtracklayer::import(paste(module_path,gtf_merged,sep=""))

# I only want exons in granges
gtf_full <- gtf_full[(elementMetadata(gtf_full)[,"type"] == "exon")]

union_of_intersect <- read.csv(paste(module_path,"masterlists/union_of_intersects_events.csv",sep=""))
union_of_intersect <- union_of_intersect[,c("Gene_name","Coordinates", "Strand")]

union_of_intersect[c('seqnames', 'ranges')] <- str_split_fixed(union_of_intersect$Coordinates, ':', 2)
union_of_intersect[c('start', 'end')] <- str_split_fixed(union_of_intersect$ranges, '-', 2)

union_of_intersect <- union_of_intersect[,c("seqnames","start","end","Strand")]

colnames(union_of_intersect) <- c("chr", "start", "end", "strand")
union_of_intersect_granges <- makeGRangesFromDataFrame(union_of_intersect,keep.extra.columns=TRUE)

# getting all exons in l that contain my high conf splicing coords
exon_hits <- gtf_full[queryHits(findOverlaps(gtf_full, union_of_intersect_granges, type="any")),] # this is to find every entry in gtf_full that contains my interested range in exon_hits

exon_hits <- as.data.frame(exon_hits)
novel_tx <- filter(exon_hits,grepl("MSTRG",transcript_id))
novel_tx <- unique(novel_tx$transcript_id)

# getting full transcript entries of the novel tx we are interested in. Gtf here only contains the novel transcripts
gtf_full <- rtracklayer::import(paste(module_path,gtf_merged,sep=""))
gtf_subset <- gtf_full[(elementMetadata(gtf_full)[,"transcript_id"] %in% novel_tx)]
rtracklayer::export(gtf_subset,paste(module_path,merged_filtered_gtf,sep=""))

# the full gtf with reference transcripts, plus the novel transcripts that we have filtered
if (species == "Mus_musculus"){
  ref_tx <- gtf_full[grepl("ENSMUST",(elementMetadata(gtf_full)[,"transcript_id"]))]
} else if (species == "Homo_sapiens"){
  ref_tx <- gtf_full[grepl("ENST",(elementMetadata(gtf_full)[,"transcript_id"]))]
}

ref_tx <- as.data.frame(ref_tx)
gtf_subset <- as.data.frame(gtf_subset)
filtered_gtf <- rbind(ref_tx,gtf_subset)
rtracklayer::export(filtered_gtf,paste(module_path,merged_fil_withRef_gtf,sep=""))

lgr$info("Filtered Stringtie file created at: %s",paste(module_path,merged_filtered_gtf,sep=""))
lgr$remove_appender("snakemake_file_log")
