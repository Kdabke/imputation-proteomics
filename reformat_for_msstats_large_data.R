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

#library(MSstats)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(bit64)

# first input name of file- txt and second is name of output file in .csv format
args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
fn_2 <- args[2]

df <- read.table(fn, stringsAsFactors = FALSE, header = T, sep = "\t")

#df <- read.csv(fn, stringsAsFactors = F)

#df <- subset(df, select = -X)

#df[, 4:23] <- lapply(4:23, function(x) as.numeric(df[[x]]))

df <- gather(df, key = "Run", value = "Intensity", -FragmentIon, -ProteinName, -PeptideSequence)

df$PrecursorCharge = str_sub(df$FragmentIon,-1)


df$random <- df$FragmentIon

df <- separate(df, col = random, into = c("a", "b", "ProductCharge","d", "e"), sep = "_")


df <- subset(df, select = -c(a,b,d,e))

df$random <- df$Run

df <- separate(df, col = random, into = c("Condition", "b", "c"), sep = "_")
df <- subset(df, select = -c(b,c))

df$BioReplicate <- df$Run

df$BioReplicate <- substring(df$BioReplicate, 5)
df$BioReplicate = substr(df$BioReplicate,1,nchar(df$BioReplicate)-2)

df$IsotopeLabelType <- "light"

df[, c(6:7)] <- lapply(c(6:7), function(x) as.integer(df[[x]]))

df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)
                                                                                       
ref_df <- read.csv("/common/dabkek/Dilution_re_analysis_swath/feat_align_all_replicates_3_norm_msstats.csv")
                                                                                       
df <- df[,colnames(ref_df)]
                                                                                       
write.csv(df, file = paste0(fn_2), row.names = FALSE)
