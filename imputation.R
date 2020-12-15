#!/bin/env Rscript
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=90G
#
#
#$ -S /hpc/apps/R/3.6.3/bin/Rscript
#$ -e /common/dabkek/original_param_imputation_msstats/re-unnorm/scratch
#$ -o /common/dabkek/original_param_imputation_msstats/re-unnorm/scratch

library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
#library(ggplot2)
library(imputeLCMD)
library(pcaMethods)
library(doParallel)
library(missForest)
set.seed(4242)
#install.packages("viridis")

#library(viridis)
#library(corrplot)
gc()
og_msstats <- read.csv("/common/dabkek/Dilution_re_analysis_swath/feat_align_all_replicates_3_unnorm_msstats.csv")

data_imputation <- as.data.frame(og_msstats[,c(1,2,4,9:10)])


data_imputation[, 1:4] <- lapply(1:4, function(x) as.character(data_imputation[[x]]))

data_imputation <- spread(data_imputation, key = Run, value = Intensity)

# removing 4 rows with all columns with missing values

test <- data_imputation
test[is.na(test)] <- 0
test3 <- test[apply(test == 0, 1, sum) >= 9, ]

#write.csv(test3, "all_missing.csv")

data_imputation_filt <- filter(data_imputation, !FragmentIon %in% test3$FragmentIon)


prot_pept <- as.data.frame(data_imputation_filt[,c(1:3)])

data_imputation_filt <- column_to_rownames(data_imputation_filt, "FragmentIon")
data_imputation_filt <- as.data.frame(data_imputation_filt[,-c(1,2)])

# zero imputation

#og_fragment_ms_0 <- data_imputation_filt

#og_fragment_ms_0[is.na(og_fragment_ms_0)] <- 0

#og_fragment_ms_0 <- cbind(prot_pept, og_fragment_ms_0)
#og_fragment_ms_0 <- subset(og_fragment_ms_0, select = -c(FragmentIon))
#og_fragment_ms_0 <- rownames_to_column(og_fragment_ms_0, "FragmentIon")

#write.table(og_fragment_ms_0, file = "og_fragment_ms_0.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# minDet

#df_minDet <- impute.MinDet(dataSet.mvs = data_imputation_filt, q = 0.0001)

#all(rownames(df_minDet) == prot_pept$FragmentIon)

#og_fragment_ms_minDet <- cbind(prot_pept, df_minDet)
#rownames(og_fragment_ms_minDet) <- NULL
#write.table(og_fragment_ms_minDet, file = "og_fragment_ms_minDet.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# minProb

#df_minProb <- impute.MinProb(dataSet.mvs = data_imputation_filt)
#all(rownames(df_minProb) == prot_pept$FragmentIon)

#og_fragment_ms_minProb <- cbind(prot_pept, df_minProb)
#rownames(og_fragment_ms_minProb) <- NULL

#write.table(og_fragment_ms_minProb, file = "og_fragment_ms_minProb.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# knn
#df_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(data_imputation_filt), K = 10))
#all(rownames(df_knn) == prot_pept$FragmentIon)

#og_fragment_ms_knn <- cbind(prot_pept, df_knn)
#rownames(og_fragment_ms_knn) <- NULL

#write.table(og_fragment_ms_knn, file = "og_fragment_ms_knn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# svd

#df_svd <- as.data.frame(imputeLCMD::impute.wrapper.SVD(data_imputation_filt, K = 8))
#all(rownames(df_svd) == prot_pept$FragmentIon)

#og_fragment_ms_svd <- cbind(prot_pept, df_svd)
#rownames(og_fragment_ms_svd) <- NULL

#write.table(og_fragment_ms_svd, file = "og_fragment_ms_svd.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# bpca 
#df_bpca <- pcaMethods::pca(as.matrix(data_imputation_filt), nPcs = 8, method = "bpca")

#data_df_bpca <- as.data.frame(completeObs(df_bpca))

#all(rownames(data_df_bpca) == prot_pept$FragmentIon)

#data_df_bpca <- cbind(prot_pept, data_df_bpca)
#rownames(data_df_bpca) <- NULL

#write.table(data_df_bpca, file = "og_fragment_ms_bpca.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# lls
#data_imputation_filt_t <- as.data.frame(t(data_imputation_filt))

#col_no <- ncol(data_imputation_filt_t)/2

#df1 <- data_imputation_filt_t[,c(1:col_no)]

#col_no2 <- col_no + 1

#df2 <- data_imputation_filt_t[,c(col_no2:ncol(data_imputation_filt_t))]


#df_lls1 <- pcaMethods::llsImpute(as.matrix(df1), k = 150, allVariables = TRUE)

#df_lls2 <- pcaMethods::llsImpute(as.matrix(df2), k = 150, allVariables = TRUE)

#result1 <- completeObs(df_lls1)
#data_df_lls1 <- as.data.frame(t(result1))

#result2 <- completeObs(df_lls2)
#data_df_lls2 <- as.data.frame(t(result2))

#all(colnames(data_df_lls1) == colnames(data_df_lls2))


#data_df_lls <- rbind(data_df_lls1, data_df_lls2)

#write.table(data_df_lls, file = "og_lls2.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# random forest

registerDoParallel(cores=8)
df_rf <- missForest(data_imputation_filt, parallelize = "forests")

data_df_rf <- df_rf$ximp
#write.table(data_df_rf, file = "data_df_rf.txt", quote = FALSE, sep = "\t", row.names = FALSE)

all(rownames(data_df_rf) == prot_pept$FragmentIon)

data_df_rf <- cbind(prot_pept, data_df_rf)
rownames(data_df_rf) <- NULL

write.table(data_df_rf, file = "og_fragment_ms_rf.txt", quote = FALSE, sep = "\t", row.names = FALSE)

