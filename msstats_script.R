#!/bin/env Rscript
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=24G
#
#
#$ -S /hpc/apps/R/4.0.2/bin/Rscript
#$ -e /common/dabkek/original_param_imputation_msstats/re-unnorm/scratch
#$ -o /common/dabkek/original_param_imputation_msstats/re-unnorm/scratch

#module load R

library(MSstats)
library(dplyr)
library(tidyr)
library(tibble)

# takes in a file name from user, a .csv file in the same working directory as first argument
# takes in name of msstats_processed file name from user as second argument
# takes in file name for protein quantification file from user as thirg argument

# expected command: qsubcommand .csv (input file name) .csv (output file name for processed file) .csv (output file name for protein quant)


args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
fn_2 <- args[2]
#fn_3 <- args[3]

df <- read.csv(fn)


df_msstats <- dataProcess(df, normalization=FALSE, MBimpute=FALSE)

#msstats_processed <- as.data.frame(df_msstats$ProcessedData, stringsAsFactors = FALSE)

sampleQuant_df <- quantification(df_msstats)

#write.csv(msstats_processed, paste0(fn_2))
write.csv(sampleQuant_df, paste0(fn_2))


