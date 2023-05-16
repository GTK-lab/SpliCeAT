# Run using command line: Rscript prep.R

library(dplyr)

########### change these ##############
setwd("/mnt/cbis/home/yongshan/whippet_snakemake")
base_path <- getwd()
fasta_file_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa.gz"
annotation_gtf_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz"
experiment_name <- "WT_IgG2A_WT_O9_CTX"
design <- read.csv("./input/design.csv")
fq_dir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/trimmedFASTQ/"
bam_dir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR/"
julia <- "/mnt/cbis/home/yongshan/julia-1.7.2/bin/julia"
whippet_bin <- "/mnt/cbis/home/yongshan/Whippet.jl/bin/"
#######################################

design$psi_gz <- paste(base_path,"/results/quantify/",design$sample,"/",design$sample,".psi.gz", sep="")

# config file

sink(paste("./config/config.yaml",sep=""))
cat(paste("samples_tsv: ",base_path,"/config/samples.tsv",sep=""))
cat("\n")
cat(paste("fasta_file: ",fasta_file_path,sep=""))
cat("\n")
cat(paste("annotation_gtf: ",annotation_gtf_path,sep=""))
cat("\n")
cat(paste("fastq_tsv: ",base_path,"/config/fastq.tsv",sep=""))
cat("\n")
cat(paste("delta_tsv: ",base_path,"/config/delta.tsv",sep=""))
cat("\n")
cat(paste("delta_input_tsv: ",base_path,"/config/delta_input.tsv",sep=""))
cat("\n")
cat(paste("bam_dir: ",bam_dir,sep=""))
cat("\n")
cat(paste("julia: ",julia,sep=""))
cat("\n")
cat(paste("whippet_bin: ",whippet_bin,sep=""))
sink()

# delta file
sink(paste("./config/delta.tsv",sep=""))
cat("experiment")
cat("\t")
cat("grp1")
cat("\t")
cat("grp2")
cat("\n")
cat(experiment_name)
cat("\t")
cat(paste(filter(design, design == "control")$psi_gz,collapse=","))
cat("\t")
cat(paste(filter(design, design == "treatment")$psi_gz,collapse=","))
sink()

# delta_input file
design$experiment <- experiment_name
write.table(design[,c("experiment","psi_gz")], file = "./config/delta_input.tsv", row.names=FALSE, sep="\t",quote=FALSE)

# fastq
fastq <- design[,c("sample","fq1","fq2")]
fastq$dir <- fq_dir
fastq$index_dir <- paste(base_path,"/results/index/",experiment_name,"/",experiment_name,".jls", sep="")
write.table(fastq, file = "./config/fastq.tsv", row.names=FALSE, sep="\t",quote=FALSE)

# samples
design$sample_alignment <- paste(design$sample, "_Aligned.sortedByCoord.out.bam", sep="")
  
sink(paste("./config/samples.tsv",sep=""))
cat("experiment")
cat("\t")
cat("cko")
cat("\t")
cat("dir")
cat("\n")
cat(experiment_name)
cat("\t")
cat(paste(filter(design, design == "treatment")$sample_alignment,collapse=" "))
cat("\t")
cat(bam_dir)
sink()