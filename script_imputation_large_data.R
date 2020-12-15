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

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
#library(ggplot2)
library(imputeLCMD)
library(pcaMethods)
library(missForest)
library(doParallel)
set.seed(4242)
#install.packages("viridis")

#library(viridis)
#library(corrplot)


df <- read.table("/common/dabkek/re-NewAcquisition/feat_align_103_norm_mapDIA.txt", header = T,
sep = "\t")
df <- subset(df, select = -c(RT))

df <- as.data.frame(df[,-c(grep("pool", colnames(df)))])
rownames(df) <- NULL
df <- column_to_rownames(df, "FragmentIon")

samples_primary_paired_tum_strm <- read.table("samples_primary_paired_tum_strm.txt", 
stringsAsFactors = F)
samples_primary_paired_tum_strm$V1 <- paste(samples_primary_paired_tum_strm$V1, "1", sep = "_")
samples_primary_tumor <- as.data.frame(samples_primary_paired_tum_strm[1:23,])
colnames(samples_primary_tumor) <- "V1"

samples_primary_stroma <- as.data.frame(samples_primary_paired_tum_strm[24:41,])
colnames(samples_primary_stroma) <- "V1"

tumor <- as.data.frame(df[,c(grep("ph", colnames(df)))])
tumor <- subset(tumor, select = c(samples_primary_tumor$V1))
stroma <- as.data.frame(df[,c(grep("strm", colnames(df)))])
stroma <- subset(stroma, select = c(samples_primary_stroma$V1))
tumor[is.na(tumor)] <- 0

tumor <- tumor[apply(tumor == 0, 1, sum) >= 22, ]

stroma[is.na(stroma)] <- 0

stroma <- stroma[apply(stroma == 0, 1, sum) >= 17, ]

df <- filter(df, !rownames(df) %in% rownames(tumor))
df <- filter(df, !rownames(df) %in% rownames(stroma))

df <- rownames_to_column(df, "FragmentIon")
prot_pept <- as.data.frame(df[,1:3])

df <- column_to_rownames(df, "FragmentIon")

df <- subset(df, select = -c(ProteinName, PeptideSequence))

df <- subset(df, select = c(samples_primary_paired_tum_strm$V1))



# bpca 
df <- as.data.frame(log2(df))
df_bpca <- pcaMethods::pca(as.matrix(df), nPcs = 10, method = "bpca")

data_df_bpca <- as.data.frame(completeObs(df_bpca))
data_df_bpca <- as.data.frame(2^data_df_bpca)

all(rownames(data_df_bpca) == prot_pept$FragmentIon)

data_df_bpca <- cbind(prot_pept, data_df_bpca)
rownames(data_df_bpca) <- NULL

write.table(data_df_bpca, file = "fragment_primary_bpca_k10_log2.txt", quote = F, sep = "\t", row.names = F)



# random forest
df <- as.data.frame(log2(df))

registerDoParallel(cores=8)
df_rf <- missForest(df, parallelize = "forests")

data_df_rf <- df_rf$ximp
data_df_rf <- as.data.frame(2^data_df_rf)

all(rownames(data_df_rf) == prot_pept$FragmentIon)

data_df_rf <- cbind(prot_pept, data_df_rf)
rownames(data_df_rf) <- NULL
write.table(data_df_rf, file = "fragment_primary_rf_log2.txt", quote = F, sep = "\t", row.names = F)




