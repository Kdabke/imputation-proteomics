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
library(stringr)

# first input name of file- txt and second is name of output file in .csv format
args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
fn_2 <- args[2]

df <- read.table(fn, stringsAsFactors = FALSE)

colnames(df) <- df[1,]
df <- as.data.frame(df[-1,])
df[, 4:12] <- lapply(4:12, function(x) as.numeric(df[[x]]))

df <- gather(df, key = "Run", value = "Intensity", -FragmentIon, -ProteinName, -PeptideSequence)

df$PrecursorCharge = str_sub(df$FragmentIon,-1)


df$random <- df$FragmentIon

df <- separate(df, col = random, into = c("a", "b", "ProductCharge","d", "e"))


df <- subset(df, select = -c(a,b,d,e))

df$random <- df$Run

df <- separate(df, col = random, into = c("Condition", "BioReplicate", "c"))
df <- subset(df, select = -c(c))
df$IsotopeLabelType <- "light"

df[, c(5:7,9)] <- lapply(c(5:7,9), function(x) as.integer(df[[x]]))

df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)
                                                                                       
ref_df <- read.csv("/common/dabkek/Dilution_re_analysis_swath/feat_align_all_replicates_3_norm_msstats.csv")
                                                                                       
df <- df[,colnames(ref_df)]
                                                                                       
write.csv(df, file = paste0(fn_2), row.names = FALSE)
