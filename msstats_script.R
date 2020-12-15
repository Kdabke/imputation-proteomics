#!/bin/env Rscript
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=90G
#
#
#$ -S /hpc/apps/R/4.0.2/bin/Rscript
#$ -e /common/dabkek/re-NewAcquisition/primary_imputation/scratch
#$ -o /common/dabkek/re-NewAcquisition/primary_imputation/scratch

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



df <- read.csv(fn)

#df <- filter(df, !Intensity < 0)

df_msstats <- dataProcess(df, normalization=FALSE, MBimpute=FALSE)

#save(df_msstats, file= "msstats_dataprocess_tum_strm_20.RData")

sampleQuant_df <- quantification(df_msstats)


write.csv(sampleQuant_df, paste0(fn_2))

