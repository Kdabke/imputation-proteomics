#!/bin/env Rscript
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=24G
#
#
#$ -S /hpc/apps/R/4.0.2/bin/Rscript
#$ -e /common/dabkek/re-NewAcquisition/prot_imputation/scratch
#$ -o /common/dabkek/re-NewAcquisition/prot_imputation/scratch




library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
#library(ComplexHeatmap)
#library(RColorBrewer)
library(ROTS)
#library(colorspace)
#library(UpSetR)

args <- commandArgs(trailingOnly = TRUE)
fn <- args[1]
fn_1 <- args[2]
#fn_2 <- args[3]
#fn_3 <- args[4]

df <- read.csv(fn)

df <- as.data.frame(df[-1,-1], stringsAsFactors = F)
df <- separate(df, col = "Protein", into = c("a","b","Protein"), sep = "\\|")
df <- subset(df, select = -c(a,b))

rownames(df) <- NULL
df <- column_to_rownames(df, "Protein")

filter <- read.table("/common/dabkek/re-NewAcquisition/prot_imputation/samples_primary_paired_tum_strm.txt")

df <- subset(df, select = c(Protein,filter$V1))

rots <- ROTS(df, groups = c(rep(0, each = 23), rep(1, each = 18)),
                           B = 1000, K = 2000, seed = 4242, paired = F, 
                           progress = TRUE, verbose = TRUE)

fdr_rots <- as.data.frame(rots$FDR)
logfc_rots <- as.data.frame(rots$logfc)
pvalue_rots <- as.data.frame(rots$pvalue)

rots <- cbind(logfc_rots, pvalue_rots, fdr_rots)

colnames(rots) <- c("log2FoldChange", "pvalue", "FDR")


rots$col_FDR <- rots$FDR <= 0.05

rots$col_pvalue <- rots$pvalue <= 0.05

write.csv(rots, paste0(fn_1))



#pdf(paste0(fn_2))


#highlight_prot <- c("PAX8", "FKBP4", "CADH1", "MSLN", "MAL2",
#                    "CAV1", "COEA1", "LUM", "EHD2", "CO6A3", "PGS2",
#                    "PGS1", "CO6A2")

#rots <- rownames_to_column(rots, X)
#rots <- separate(rots, col = "X", into = "Protein", sep = "_")

#pdf(paste0(fn_2))

#ggplot(data = rots, aes(x = log2FoldChange, y = -1*log10(pvalue))) + geom_point(size = 3, color = "grey") +
#    geom_point(data = filter(rots, log2FoldChange >= 0.5) %>% filter(., FDR <= 0.05), color = "#e5f5f9", size = 4) +
#    geom_point(data = filter(rots, log2FoldChange <= -0.5) %>% filter(., FDR <= 0.05), color = "#fee8c8", size = 4) +
#    geom_point(data = filter(rots, Protein %in% highlight_prot), color = "blue", size = 3) +
#    geom_text(data = filter(rots, Protein %in% highlight_prot), aes(label = Protein), nudge_x = -0.10, nudge_y = 0.05, size = 3) +
#    labs(y = "-log10 pvalue") +
#    theme_bw() +
#    theme(panel.grid = element_blank(),
#          text = element_text(size = 20)) +
#    labs(title = fn_3) 


#dev.off()


