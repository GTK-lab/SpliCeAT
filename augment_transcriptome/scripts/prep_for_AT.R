#### CHANGE THIS ####
config_file_path <- "/mnt/cbis/home/yongshan/SpliCeAT/augment_transcriptome/config/config.yaml"
####################

library(lgr)

# libraries
lgr$info("Loading libraries...")

suppressMessages({
library(yaml)
library(biomaRt)
library(stringr)
library(data.table)
library(VennDetail)
suppressMessages(library(dplyr))
library(GenomicRanges)
})

lgr$info("Done.")

# load in config containing file paths and other params
lgr$info("Reading config.yaml...")

config <- read_yaml(config_file_path)

lgr$info("Done.")

###################
# MASTERLISTS PREP
###################

lgr$info("MASTERLISTS PREP")

# Prepares the Majiq, Whippet and Leafcutter splicing results into a masterlist for comparison purposes

# MAJIQ
lgr$info("Getting Ensembl annotations...")

ensembl <- useEnsembl(biomart = "genes",
                      dataset = "mmusculus_gene_ensembl",
                      version = config$mouse_ensembl_version)

annotations <- getBM(attributes = c('ensembl_gene_id',"external_gene_name",'description', 'chromosome_name',
                                    'start_position', 'end_position', 'strand'), mart = ensembl)

lgr$info("Done.")

lgr$info("Processing Majiq LSV file...")

majiq_lsv <- read.table(config$majiq_lsv_file_path)
majiq_lsv[c('Ensembl_ID', 'type',"Coord")] <- str_split_fixed(majiq_lsv$V1, ':', 3)
majiq_lsv <- majiq_lsv[,c(2,4)]
gene_ids <- str_replace(majiq_lsv$Ensembl_ID,
                        pattern = ".[0-9]+$",
                        replacement = "")
majiq_lsv <- cbind(gene_ids,majiq_lsv)
majiq_lsv_final <- left_join(majiq_lsv,
                             annotations,
                             by = c("gene_ids"="ensembl_gene_id"))

majiq_lsv_final <- majiq_lsv_final[,c(1:4)]

dir.create(paste(config$BASE_PATH,"masterlists",sep=""))
write.csv(majiq_lsv_final,paste(config$BASE_PATH,"masterlists/majiq_masterlist.csv",sep=""),row.names = FALSE)

lgr$info("Majiq masterlist created at: %s",paste(config$BASE_PATH,"masterlists/majiq_masterlist.csv",sep=""))

# LEAFCUTTER
lgr$info("Processing Leafcutter input files...")
leafcutter_cluster_sig <- read.table(config$leafcutter_cluster_sig_file_path,sep="\t",header=TRUE)
leafcutter_cluster_sig <- leafcutter_cluster_sig[,c(1,6,7)]
leafcutter_effect_sizes <- read.table(config$leafcutter_effect_size_file_path,sep="\t",header=TRUE)
leafcutter_effect_sizes[c("chr","start","stop","cluster")] <- str_split_fixed(leafcutter_effect_sizes$intron, ':', 4)
leafcutter_effect_sizes$cluster_joined <- paste(leafcutter_effect_sizes$chr,":",leafcutter_effect_sizes$cluster,sep="")
leafcutter_effect_sizes <- leafcutter_effect_sizes[,c(1,5,10)]

leafcutter_final <- left_join(leafcutter_effect_sizes,
                              leafcutter_cluster_sig,
                              by = c("cluster_joined"="cluster"))

leafcutter_final_filtered <- filter(leafcutter_final, p.adjust <0.05 & abs(deltapsi) >= 0.2)
write.csv(leafcutter_final_filtered,paste(config$BASE_PATH,"masterlists/leafcutter_masterlist.csv",sep=""),row.names = FALSE)

lgr$info("Leafcutter masterlist created at: %s",paste(config$BASE_PATH,"masterlists/leafcutter_masterlist.csv",sep=""))

# WHIPPET
lgr$info("Processing Whippet input file...")

whippet_diff <- read.table(config$whippet_diff_file_path, sep = '\t', header = TRUE, row.names=NULL)
colnames(whippet_diff) <- c("Gene", "Node", "Coord", "Strand", "Type", "Psi_A", "Psi_B", "DeltaPsi", "Probability", "Complexity", "Entropy")
whippet_diff <- whippet_diff[,1:11]
whippet_diff <- whippet_diff[,c(1,3,5,8,9)]
gene_ids <- str_replace(whippet_diff$Gene,
                        pattern = ".[0-9]+$",
                        replacement = "")
whippet_diff <- cbind(gene_ids,whippet_diff)
group <- as.data.table(whippet_diff)
whippet_diff_highest_prob <- group[group[, .I[Probability == max(Probability)], by=Gene]$V1]
whippet_diff_highest_prob <- as.data.frame(whippet_diff_highest_prob)
whippet_diff_highest_prob[c('chr', 'start')] <- str_split_fixed(whippet_diff_highest_prob$Coord, ':', 2)
whippet_diff_highest_prob[c('start', 'stop')] <- str_split_fixed(whippet_diff_highest_prob$start, '-', 2)

whippet_diff_final <- left_join(
  whippet_diff_highest_prob,
  annotations,
  by = c("gene_ids"="ensembl_gene_id")
)
colnames(whippet_diff_final) <- c("Ensembl_ID","Ensembl_ID_version",
                                        "CE_coord","Type",
                                        "DeltaPSI","Probability",
                                        "CE_chr", "CE_start",
                                        "CE_stop", "Gene_name",
                                        "Description","Gene_chr",
                                        "Gene_start","Gene_stop","Gene_strand")

# replacing -1 and 1 with - and +
whippet_diff_final$Gene_strand[whippet_diff_final$Gene_strand=="-1"]<-"-"
whippet_diff_final$Gene_strand[whippet_diff_final$Gene_strand=="1"]<-"+"

write.csv(whippet_diff_final, paste(config$BASE_PATH,"masterlists/whippet_masterlist.csv",sep=""),row.names = FALSE)

lgr$info("Whippet masterlist created at: %s",paste(config$BASE_PATH,"masterlists/whippet_masterlist.csv",sep=""))

# Getting intersections of majiq, whippet, leafcutter
lgr$info("Getting intersections of Majiq, Whippet, Leafcutter...")

filter_genes <- function(whippet_csv, majiq_csv, leafcutter_csv){
  whippet_df <- read.csv(whippet_csv,header=TRUE,check.names=FALSE)
  majiq_df <- read.csv(majiq_csv,header=TRUE,check.names=FALSE)
  leafcutter_df <- read.csv(leafcutter_csv,header=TRUE,check.names=FALSE)
  
  whippet_df_filtered <- filter(whippet_df, abs(DeltaPSI) >= 0.2 & Probability >= 0.95)
  whippet_df_filtered <- whippet_df_filtered[!(is.na(whippet_df_filtered$Gene_name) | whippet_df_filtered$Gene_name==""), ]
  whippet_genes <- sort(unique(na.omit(whippet_df_filtered[["Gene_name"]])))
  
  majiq_df <- majiq_df[!(is.na(majiq_df$external_gene_name) | majiq_df$external_gene_name==""), ]
  majiq_genes <- sort(unique(na.omit(majiq_df[["external_gene_name"]])))
  
  leafcutter_df <- leafcutter_df[!(is.na(leafcutter_df$genes) | leafcutter_df$genes==""), ]
  leafcutter_genes <- sort(unique(na.omit(leafcutter_df[["genes"]])))
  
  result <- list("majiq_genes" = majiq_genes, "whippet_genes" = whippet_genes, "leafcutter_genes" = leafcutter_genes)
}

filtered_genes <- filter_genes(paste(config$BASE_PATH,"masterlists/whippet_masterlist.csv",sep=""),
                               paste(config$BASE_PATH,"masterlists/majiq_masterlist.csv",sep=""),
                               paste(config$BASE_PATH,"masterlists/leafcutter_masterlist.csv",sep=""))

master_whippet_df <- read.csv(paste(config$BASE_PATH,"masterlists/whippet_masterlist.csv",sep=""),header=TRUE,check.names=FALSE)

master_whippet_df$Description <- sub(" \\[.*", "", master_whippet_df$Description)

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
Leafcutter <- as.data.frame(getSet(ven, subset = c("Leafcutter"))$Detail)
colnames(Leafcutter) <- "Gene_name"
Whippet_MAJIQ <- as.data.frame(getSet(ven, subset = c("Whippet_MAJIQ"))$Detail)
colnames(Whippet_MAJIQ) <- "Gene_name"
MAJIQ <- as.data.frame(getSet(ven, subset = c("MAJIQ"))$Detail)
colnames(MAJIQ) <- "Gene_name"
Whippet <- as.data.frame(getSet(ven, subset = c("Whippet"))$Detail)
colnames(Whippet) <- "Gene_name"

# only want union of intersect

union <- rbind(shared, MAJIQ_Leafcutter,Whippet_Leafcutter,Whippet_MAJIQ)
union <- as.data.frame(union[order(union$Gene_name),])
colnames(union) <- "Gene_name"

df_to_display <- merge(union, master_whippet_df, by = 'Gene_name')
df_to_display <- df_to_display[,c(1,11,4,5,6,7,15)]
colnames(df_to_display) <- c("Gene", "Description",
                             "Coordinates", "Type",
                             "DeltaPSI", "Probability",
                             "Strand")

write.csv(df_to_display, paste(config$BASE_PATH,"masterlists/union_of_intersects_events.csv",sep=""),row.names = FALSE)

lgr$info("Intersections file created at: %s",paste(config$BASE_PATH,"masterlists/union_of_intersects_events.csv",sep=""))

#############################
# FILTERING FOR SPLICE EVENTS
#############################

lgr$info("FILTERING FOR SPLICE EVENTS")
lgr$info("Processing merged Stringtie GTF file...")

gtf_full <- rtracklayer::import(paste(config$BASE_PATH,"results/merged_assembly/merged_stringtie_assembly.gtf",sep=""))

# I only want exons in granges
gtf_full <- gtf_full[(elementMetadata(gtf_full)[,"type"] == "exon")]

union_of_intersect <- read.csv(paste(config$BASE_PATH,"masterlists/union_of_intersects_events.csv",sep=""))
union_of_intersect <- union_of_intersect[,c("Coordinates", "Strand")]
union_of_intersect[c('seqnames', 'ranges')] <- str_split_fixed(union_of_intersect$Coordinates, ':', 2)
union_of_intersect[c('start', 'end')] <- str_split_fixed(union_of_intersect$ranges, '-', 2)
union_of_intersect <- union_of_intersect[,c("seqnames","start","end","Strand")]
colnames(union_of_intersect) <- c("chr", "start", "end", "strand")
union_of_intersect_granges <- makeGRangesFromDataFrame(union_of_intersect,keep.extra.columns=TRUE)

# getting all exons in gtf_full that contain my high conf splicing coords
exon_hits <- gtf_full[queryHits(findOverlaps(gtf_full, union_of_intersect_granges, type="any")),] # this is to find every entry in gtf_full that contains my interested range in exon_hits
exon_hits <- as.data.frame(exon_hits)
novel_tx <- filter(exon_hits,grepl("MSTRG",transcript_id))
novel_tx <- unique(novel_tx$transcript_id)

# getting full transcript entries of the novel tx we are interested in. Gtf here only contains the novel transcripts
gtf_full <- rtracklayer::import(paste(config$BASE_PATH,"results/merged_assembly/merged_stringtie_assembly.gtf",sep=""))
gtf_subset <- gtf_full[(elementMetadata(gtf_full)[,"transcript_id"] %in% novel_tx)]
rtracklayer::export(gtf_subset,paste(config$BASE_PATH,"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf",sep=""))

# the full gtf with reference transcripts, plus the novel transcripts that we have filtered
ref_tx <- gtf_full[grepl("ENSMUST",(elementMetadata(gtf_full)[,"transcript_id"]))]
ref_tx <- as.data.frame(ref_tx)
gtf_subset <- as.data.frame(gtf_subset)
filtered_gtf <- rbind(ref_tx,gtf_subset)
rtracklayer::export(filtered_gtf,paste(config$BASE_PATH,"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered_with_reference.gtf",sep=""))

lgr$info("Filtered Stringtie file created at: %s",paste(config$BASE_PATH,"results/merged_assembly/merged_stringtie_assembly_novel_exon_filtered.gtf",sep=""))

