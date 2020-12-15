#!/bin/env Rscript

#module load R

library(dplyr)
library(tidyr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
fn_2 <- args[2]
fn_3 <- args[3]

df <- read.csv(fn, stringsAsFactors = FALSE)

df <- as.data.frame(df[-1,-1])
#rownames(df) <- NULL

#df <- column_to_rownames(df, "Protein")

df[is.na(df)] <- 0

#df <- as.data.frame(log2(df))

#df <- rownames_to_column(df, "Protein")


df$high_mean <- rowMeans(df[,2:4])
df$low_mean <- rowMeans(df[,5:7])
df$med_mean <- rowMeans(df[,8:10])

df_fldchnge <- as.data.frame(df[,c(1,11:13)])

df_fldchnge$high_med <- df_fldchnge$high_mean - df_fldchnge$med_mean
df_fldchnge$med_low <- df_fldchnge$med_mean - df_fldchnge$low_mean
df_fldchnge$high_low <- df_fldchnge$high_mean - df_fldchnge$low_mean

df_fldchnge <- as.data.frame(df_fldchnge[,-c(2:4)])

df_fldchnge <- gather(df_fldchnge, key = "load", value = "foldchange", -Protein)

df_fldchnge$imputation <- fn_2

write.table(df_fldchnge, file = paste0(fn_3), quote = FALSE, sep = "\t", row.names = FALSE)
