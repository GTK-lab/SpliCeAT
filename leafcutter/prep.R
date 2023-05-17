# Run using command line: Rscript prep.R

########### change these ##############
setwd("/mnt/cbis/home/yongshan/leafcutter_snakemake")
base_path <- getwd()
experiment_name <- "WT_IgG2A_WT_O9_CTX"
design <- read.csv("./input/design.csv")
regtools_strand <- "1"
gene_annotation <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz"
bam_dir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR/"
leafcutter_dir <- "/mnt/cbis/home/yongshan/leafcutter/"
#######################################

# XXX_groups_file
write.table(design, file = paste(base_path,"/config/",experiment_name,"_groups_file.txt",sep=""), 
            row.names=FALSE, sep="\t",quote=FALSE,col.names=FALSE)

# XXX_junc
design$junc <- paste(base_path, "/results/", design$sample, ".junc", sep="")
write.table(design$junc, file = paste(base_path,"/config/",experiment_name,"_junc.txt",sep=""), 
            row.names=FALSE, sep="\t",quote=FALSE,col.names=FALSE)

# config
sink(paste("./config/config.yaml",sep=""))
cat(paste("BASE_PATH: ",base_path,"/",sep=""))
cat("\n")
cat(paste("samples_tsv: ",base_path,"/config/samples.tsv",sep=""))
cat("\n")
cat(paste("STRAND: ",regtools_strand,sep=""))
cat("\n")
cat(paste("juncs_file: ",base_path,"/config/juncs_file.tsv",sep=""))
cat("\n")
cat(paste("output_juncs: ",base_path,"/config/output_junc.tsv",sep=""))
cat("\n")
cat(paste("gtf: ",gene_annotation,sep=""))
cat("\n")
cat(paste("leafcutter_dir: ",leafcutter_dir,sep=""))
cat("\n")
cat(paste("bam_dir: ",bam_dir,sep=""))
sink()

# juncs_file
sink(paste("./config/juncs_file.tsv",sep=""))
cat("experiment")
cat("\t")
cat("juncs_file")
cat("\n")
cat(experiment_name)
cat("\t")
cat(paste(experiment_name,"_junc.txt",sep=""))
sink()

# output_junc
design$experiment <- experiment_name
design$junction_files <- paste(design$sample,".junc",sep="")
write.table(design[,c("experiment","junction_files")],
            file=paste(base_path,"/config/output_junc.tsv",sep=""),
            row.names=FALSE, sep="\t", quote=FALSE)

# samples
design$folder <- experiment_name
design$bam <- paste(design$sample,"_Aligned.sortedByCoord.out.bam",sep="")
write.table(design[,c("sample","folder","bam")],
            file=paste(base_path,"/config/samples.tsv",sep=""),
            row.names=FALSE, sep="\t",quote=FALSE)
