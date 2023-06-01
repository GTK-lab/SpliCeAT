#!/usr/bin/env Rscript

library(dplyr)

############################################### change these ##################################################
# differential splicing detection directory
setwd("/mnt/cbis/home/yongshan/SpliCeAT/ds_detection") 

# experimental design csv file - dont need to change if design.csv is in input directory
design <- read.csv("./input/design.csv") 

# directory containing your bam alignment files
bam_dir <- "/mnt/gtklab01/linglab/tdp43/STAR/tdp43_nestin_ctx_e14/" 

 # directory containing your fastq files
fq_dir <- "/mnt/gtklab01/linglab/tdp43/fastq/"

# set your own experiment name
experiment_name <- "tdp43_nestin_ctx_e14" 

# path of annotation gff3 file
gff3_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gff3" 

# path of annotation gtf.gz file
annotation_gtf_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM29.primary_assembly.annotation.gtf.gz" 

# path of genome fa.gz file
fasta_file_path <- "/mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa.gz" 

# path of julia command
julia <- "/mnt/cbis/home/yongshan/julia-1.7.2/bin/julia" 

# directory of whippet scripts
whippet_bin <- "/mnt/cbis/home/yongshan/Whippet.jl/bin/" 

# Note that: 0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR
regtools_strand <- "1" 

# directory of leafcutter installation
leafcutter_dir <- "/mnt/cbis/home/yongshan/leafcutter/" 
###############################################################################################################

base_path <- getwd()

# config
sink(paste("./config/config.yaml",sep=""))
cat(paste("confs: ",base_path,"/config/confs.tsv",sep=""))
cat("\n")
cat(paste("delta_psi_samples: ",base_path,"/config/delta_psi_samples.tsv",sep=""))
cat("\n")
cat(paste("BASE_PATH: ",base_path,sep=""))
cat("\n")
cat(paste("gff3: ",gff3_path,sep=""))
cat("\n")
cat(paste("build_out: ",base_path,"/results/majiq_build/{experiment}",sep=""))
cat("\n")
cat(paste("delta_out: ",base_path,"/results/majiq_delta_psi/{experiment}",sep=""))
cat("\n")
cat(paste("experiment_sample_names: ",base_path,"/config/experiment_sample_names.tsv",sep=""))
cat("\n")
cat(paste("samples_whippet_tsv: ",base_path,"/config/samples_whippet.tsv",sep=""))
cat("\n")
cat(paste("samples_leafcutter_tsv: ",base_path,"/config/samples_leafcutter.tsv",sep=""))
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
cat("\n")
cat(paste("STRAND: ",regtools_strand,sep=""))
cat("\n")
cat(paste("juncs_file: ",base_path,"/config/juncs_file.tsv",sep=""))
cat("\n")
cat(paste("output_juncs: ",base_path,"/config/output_junc.tsv",sep=""))
cat("\n")
cat(paste("leafcutter_dir: ",leafcutter_dir,sep=""))
sink()

# majiq files

design$sample_alignment <- paste(design$sample, "_Aligned.sortedByCoord.out", sep="")

# conf file
sink(paste("./conf/",experiment_name,"_conf.txt",sep=""))
cat("[info]")
cat("\n")
cat(paste("bamdirs=",bam_dir,sep=""))
cat("\n")
cat(paste("genome=",unique(design$genome),sep=""))
cat("\n")
cat(paste("strandness=",unique(design$strand),sep=""))
cat("\n")
cat("[experiments]")
cat("\n")
cat(paste(experiment_name,"=",paste(design$sample_alignment,collapse=","),sep=""))
sink()

# confs file
sink(paste("./config/","confs.tsv",sep=""))
cat("experiment")
cat("\t")
cat("conf_path")
cat("\n")
cat(experiment_name)
cat("\t")
cat(paste(base_path,"/conf/",experiment_name,"_conf.txt",sep=""))
sink()

# delta_psi_samples file
sink(paste("./config/delta_psi_samples.tsv",sep=""))
cat("experiment")
cat("\t")
cat("grp1")
cat("\t")
cat("grp2")
cat("\t")
cat("dir")
cat("\n")
cat(experiment_name)
cat("\t")
cat(paste(filter(design, design == "control")$sample_alignment,collapse=" "))
cat("\t")
cat(paste(filter(design, design == "treatment")$sample_alignment,collapse=" "))
cat("\t")
cat(paste(base_path,"/results/majiq_build/",experiment_name,"/",sep=""))
sink()

# experiment_sample_names file
design_experiment_sample_names <- design[,c("sample","sample_alignment")]
colnames(design_experiment_sample_names) <- c("experiment","sample")
write.table(design_experiment_sample_names, file = "./config/experiment_sample_names.tsv", row.names=FALSE, sep="\t",quote=FALSE)

# whippet files

design$psi_gz <- paste(base_path,"/results/quantify/",design$sample,"/",design$sample,".psi.gz", sep="")

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

# samples_whippet.tsv
design$sample_alignment <- paste(design$sample, "_Aligned.sortedByCoord.out.bam", sep="")
  
sink(paste("./config/samples_whippet.tsv",sep=""))
cat("experiment")
cat("\t")
cat("treatment")
cat("\t")
cat("dir")
cat("\n")
cat(experiment_name)
cat("\t")
cat(paste(filter(design, design == "treatment")$sample_alignment,collapse=" "))
cat("\t")
cat(bam_dir)
sink()

# leafcutter files

# XXX_groups_file
write.table(design, file = paste(base_path,"/config/",experiment_name,"_groups_file.txt",sep=""), 
            row.names=FALSE, sep="\t",quote=FALSE,col.names=FALSE)
			
# XXX_junc
design$junc <- paste(base_path, "/results/", design$sample, ".junc", sep="")
write.table(design$junc, file = paste(base_path,"/config/",experiment_name,"_junc.txt",sep=""), 
            row.names=FALSE, sep="\t",quote=FALSE,col.names=FALSE)

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
			
# samples_leafcutter.tsv
design$folder <- experiment_name
design$bam <- paste(design$sample,"_Aligned.sortedByCoord.out.bam",sep="")
write.table(design[,c("sample","folder","bam")],
            file=paste(base_path,"/config/samples_leafcutter.tsv",sep=""),
            row.names=FALSE, sep="\t",quote=FALSE)
