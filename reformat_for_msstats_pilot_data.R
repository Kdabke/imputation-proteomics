#!/bin/env Rscript
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=24G
#
#
#$ -S /hpc/apps/R/4.0.2/bin/Rscript
#$ -e /common/dabkek/rf_pocroc/impute_log2/scratch
#$ -o /common/dabkek/rf_pocroc/impute_log2/scratch

#module load R

#library(MSstats)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(bit64)

# give two arguments as input - first the input file name, second the name of the output file
args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
fn_2 <- args[2]

df <- read.table(fn, stringsAsFactors = FALSE, header = T)

df[, 4:17] <- lapply(4:17, function(x) as.numeric(df[[x]]))

df <- gather(df, key = "Run", value = "Intensity", -FragmentIon, -ProteinName, -PeptideSequence)

df$PrecursorCharge = str_sub(df$FragmentIon,-1)


df$random <- df$FragmentIon

df <- separate(df, col = random, into = c("a", "b", "ProductCharge","d", "e"))


df <- subset(df, select = -c(a,b,d,e))

df$random <- df$Run

df <- separate(df, col = random, into = c("c", "Condition", "BioReplicate"))
#df$BioReplicate <- paste(df$BioReplicate, df$c, sep = "_")

df <- subset(df, select = -c(c))
df$IsotopeLabelType <- "light"

df[, c(5:7,9)] <- lapply(c(5:7,9), function(x) as.integer64(df[[x]]))

df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)
                                                                                       
ref_df <- read.csv("/common/dabkek/Dilution_re_analysis_swath/feat_align_all_replicates_3_unnorm_msstats.csv")
                                                                                       
df <- df[,colnames(ref_df)]
                                                                                       
write.csv(df, file = paste0(fn_2), row.names = FALSE)
