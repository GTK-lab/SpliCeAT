# Run using command line: Rscript prep.R

library(dplyr)

########### change these ##############
setwd("/mnt/cbis/home/yongshan/majiq_snakemake")
design <- read.csv("./input/design.csv")
bamdir <- "/mnt/gtklab01/linglab/external_datasets/tdp43_Q331K_rescue_rubychen/STAR"
experiment_name <- "WT_IgG2A_WT_O9_CTX"
gff3 <- "gencode.vM29.primary_assembly.annotation.gff3"
base_path <- getwd()
#######################################

design$sample_alignment <- paste(design$sample, "_Aligned.sortedByCoord.out", sep="")

# conf file
sink(paste("./conf/",experiment_name,"_conf.txt",sep=""))
cat("[info]")
cat("\n")
cat(paste("bamdirs=",bamdir,sep=""))
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
design_experiment_sample_names <- design[,c(1,5)]
colnames(design_experiment_sample_names) <- c("experiment","sample")
write.table(design_experiment_sample_names, file = "./config/experiment_sample_names.tsv", row.names=FALSE, sep="\t",quote=FALSE)

# config file
sink(paste("./config/config.yaml",sep=""))
cat(paste("confs: ",base_path,"/config/confs.tsv",sep=""))
cat("\n")
cat(paste("delta_psi_samples: ",base_path,"/config/delta_psi_samples.tsv",sep=""))
cat("\n")
cat(paste("BASE_PATH: ",base_path,sep=""))
cat("\n")
cat(paste("gtf: ",base_path,"/input/",gff3,sep=""))
cat("\n")
cat(paste("build_out: ",base_path,"/results/majiq_build/{experiment}",sep=""))
cat("\n")
cat(paste("delta_out: ",base_path,"/results/majiq_delta_psi/{experiment}",sep=""))
cat("\n")
cat(paste("experiment_sample_names: ",base_path,"/config/experiment_sample_names.tsv",sep=""))
sink()
