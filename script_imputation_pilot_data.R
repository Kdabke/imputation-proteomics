#!/bin/env Rscript
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=90G
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
#library(doParallel)
#library(missForest)
library(imputeLCMD)
library(pcaMethods)
set.seed(4242)

gc()

og_msstats <- read.csv("/common/dabkek/rf_pocroc/featurealignment/feat_align_pocroc_replicate_norm_msstats.csv")

data_imputation <- as.data.frame(og_msstats[,c(1,2,4,9:10)])


data_imputation[, 1:4] <- lapply(1:4, function(x) as.character(data_imputation[[x]]))

data_imputation <- spread(data_imputation, key = Run, value = Intensity)

data_imputation <- subset(data_imputation, select = -c(pool_4_1, pool_5_1, pool_6_1))


# removing 1217 rows with all columns with missing values

test <- data_imputation
test[is.na(test)] <- 0
test3 <- test[apply(test == 0, 1, sum) >= 14, ]
write.table(test3, file = "all_missing_14_pocroc.txt", quote = FALSE, sep = "\t", row.names = FALSE)


data_imputation_filt <- filter(data_imputation, !FragmentIon %in% test3$FragmentIon)


prot_pept <- as.data.frame(data_imputation_filt[,c(1:3)])

data_imputation_filt <- column_to_rownames(data_imputation_filt, "FragmentIon")
data_imputation_filt <- as.data.frame(data_imputation_filt[,-c(1,2)])


# zero imputation 

#zero_df <- data_imputation_filt

#zero_df[is.na(zero_df)] <- 0

#all(rownames(zero_df) == prot_pept$FragmentIon)

#zero_df <- cbind(prot_pept, zero_df)
#rownames(zero_df) <- NULL
#write.table(fragment_zero_df, file = "fragment_ms_zero_norm.txt", quote = FALSE, sep = "\t", row.names = FALSE)



# minDet imputation
#df_minDet <- impute.MinDet(dataSet.mvs = data_imputation_filt, q = 0.00001)

#all(rownames(df_minDet) == prot_pept$FragmentIon)

#df_minDet <- cbind(prot_pept, df_minDet)
#rownames(og_fragment_ms_minDet) <- NULL
#write.table(og_fragment_ms_minDet, file = "pocroc_replicate_fragment_ms_minDet_norm.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# knn imputation

#df_knn <- as.data.frame(impute.wrapper.KNN(as.matrix(data_imputation_filt), K = 10))
#all(rownames(df_knn) == prot_pept$FragmentIon)

#og_fragment_ms_knn <- cbind(prot_pept, df_knn)
#rownames(og_fragment_ms_knn) <- NULL

#write.table(og_fragment_ms_knn, file = "pocroc_replicate_fragment_ms_knn_norm.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# minProb

#df_minProb <- impute.MinProb(dataSet.mvs = data_imputation_filt, tune.sigma = 0.001)
#all(rownames(df_minProb) == prot_pept$FragmentIon)

#df_minProb <- cbind(prot_pept, df_minProb)
#rownames(df_minProb) <- NULL

#write.table(df_minProb, file = "fragment_minProb.txt", quote = FALSE, sep = "\t", row.names = F)

# svd

#data_imputation_filt <- log2(data_imputation_filt)
#df_svd <- as.data.frame(imputeLCMD::impute.wrapper.SVD(data_imputation_filt, K = 10))
#df_svd <- 2^df_svd

#all(rownames(df_svd) == prot_pept$FragmentIon)

#df_svd <- cbind(prot_pept, df_svd)
#rownames(df_svd) <- NULL


#write.table(df_svd, file = "fragment_svd.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# bpca
data_imputation_filt <- log2(data_imputation_filt)
df_bpca <- pcaMethods::pca(as.matrix(data_imputation_filt), nPcs = 10, method = "bpca")

data_df_bpca <- as.data.frame(completeObs(df_bpca))
data_df_bpca <- 2^data_df_bpca

all(rownames(data_df_bpca) == prot_pept$FragmentIon)

data_df_bpca <- cbind(prot_pept, data_df_bpca)
rownames(data_df_bpca) <- NULL

write.table(data_df_bpca, file = "fragment_bpca.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# lls imputation
# splitting data set into two parts by peptide sequence
uniq_peptide <- unique(prot_pept$PeptideSequence)
fragm1 <- filter(prot_pept, PeptideSequence %in% uniq_peptide[1:5911])
fragm2 <- filter(prot_pept, PeptideSequence %in% uniq_peptide[5912:11822])
fragm3 <- filter(prot_pept, PeptideSequence %in% uniq_peptide[11823:17733])
fragm4 <- filter(prot_pept, PeptideSequence %in% uniq_peptide[17733:23644])

df1 <- filter(data_imputation_filt, rownames(data_imputation_filt) %in% c(fragm1$FragmentIon))


df1_t <- as.matrix(t(log2(df1)))

df1_lls <- pcaMethods::llsImpute(df1_t, k = 150, allVariables=TRUE, verbose = TRUE)


result1 <- completeObs(df1_lls)
data_df_lls1 <- as.data.frame(t(result1))

all(rownames(data_df_lls1) == fragm1$FragmentIon)

data_df_lls1 <- cbind(fragm1, data_df_lls1)
rownames(data_df_lls1) <- NULL

write.table(data_df_lls1, file = "data_df_lls1.txt", quote = FALSE, sep = "\t", row.names = FALSE)


df1 <- filter(data_imputation_filt, rownames(data_imputation_filt) %in% c(fragm2$FragmentIon))


df1_t <- as.matrix(t(log2(df1)))

df1_lls <- pcaMethods::llsImpute(df1_t, k = 150, allVariables=TRUE, verbose = TRUE)


result1 <- completeObs(df1_lls)
data_df_lls1 <- as.data.frame(t(result1))

all(rownames(data_df_lls1) == fragm2$FragmentIon)

data_df_lls1 <- cbind(fragm2, data_df_lls1)
rownames(data_df_lls1) <- NULL

write.table(data_df_lls1, file = "data_df_lls2.txt", quote = FALSE, sep = "\t", row.names = FALSE)

df1 <- filter(data_imputation_filt, rownames(data_imputation_filt) %in% c(fragm1[c(1:18),3], fragm3$FragmentIon))


df1_t <- as.matrix(t(log2(df1)))

df1_lls <- pcaMethods::llsImpute(df1_t, k = 150, allVariables=TRUE, verbose = TRUE)


result1 <- completeObs(df1_lls)
data_df_lls1 <- as.data.frame(t(result1))

data_df_lls1 <- filter(data_df_lls1, rownames(data_df_lls1) %in% fragm3$FragmentIon)

all(rownames(data_df_lls1) == fragm3$FragmentIon)

data_df_lls1 <- cbind(fragm3, data_df_lls1)
rownames(data_df_lls1) <- NULL

write.table(data_df_lls1, file = "data_df_lls3.txt", quote = FALSE, sep = "\t", row.names = FALSE)

df1 <- filter(data_imputation_filt, rownames(data_imputation_filt) %in% c(fragm1[c(1:18),3], fragm4$FragmentIon))


df1_t <- as.matrix(t(log2(df1)))

df1_lls <- pcaMethods::llsImpute(df1_t, k = 150, allVariables=TRUE, verbose = TRUE)


result1 <- completeObs(df1_lls)
data_df_lls1 <- as.data.frame(t(result1))

data_df_lls1 <- filter(data_df_lls1, rownames(data_df_lls1) %in% fragm4$FragmentIon)

all(rownames(data_df_lls1) == fragm4$FragmentIon)

data_df_lls1 <- cbind(fragm4, data_df_lls1)
rownames(data_df_lls1) <- NULL

write.table(data_df_lls1, file = "data_df_lls4.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# RF
registerDoParallel(cores=8)
df_rf <- missForest(log2(data_imputation_filt), parallelize = "forests")

data_df_rf <- df_rf$ximp
data_df_rf <- 2^data_df_rf
#write.table(data_df_rf, file = "data_df_rf.txt", quote = FALSE, sep = "\t", row.names = FALSE)

all(rownames(data_df_rf) == prot_pept$FragmentIon)

data_df_rf <- cbind(prot_pept, data_df_rf)
rownames(data_df_rf) <- NULL

write.table(data_df_rf, file = "fragment_rf_log2.txt", quote = FALSE, sep = "\t", row.names = FALSE)






















