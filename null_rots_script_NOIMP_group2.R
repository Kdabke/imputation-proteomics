#!/bin/env Rscript
#
#$ -cwd
#$ -pe smp 6
#$ -l mem_free=24G
#
#
#$ -S /hpc/apps/R/4.0.2/bin/Rscript
#$ -e /common/dabkek/re-NewAcquisition/primary_imputation/mapDIA/mapDIA_noimp_results_SDF2_minobs_1/scratch
#$ -o /common/dabkek/re-NewAcquisition/primary_imputation/mapDIA/mapDIA_noimp_results_SDF2_minobs_1/scratch




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

df <- read.csv(fn, stringsAsFactors = F)

df <- as.data.frame(df[-1,-1], stringsAsFactors = F)
df <- separate(df, col = "Protein", into = c("a","b","Protein"), sep = "\\|")
df <- subset(df, select = -c(a,b))

rownames(df) <- NULL
df <- column_to_rownames(df, "Protein")

set.seed(4242)
x <- as.data.frame(sample(c(0,1), 41, replace = TRUE, prob = c(0.3, 0.7)))
colnames(x) <- "group"
x$colname <- colnames(df)
x <- x[order(x$group),]

df <- as.data.frame(df[,x$colname])
grp1_data <- as.data.frame(df[,1:table(x$group)[[1]]])
grp2_data <- as.data.frame(df[,table(x$group)[[1]]+1:table(x$group)[[2]]])

grp1_data[is.na(grp1_data)] <- 0
grp2_data[is.na(grp2_data)] <- 0

grp1_data <- grp1_data[apply(grp1_data == 0, 1, sum) >= (table(x$group)[[1]]-1), ]
grp2_data <- grp2_data[apply(grp2_data == 0, 1, sum) >= (table(x$group)[[2]]-1), ]

df <- filter(df, !rownames(df) %in% rownames(grp1_data))
df <- filter(df, !rownames(df) %in% rownames(grp2_data))

rots <- ROTS(df, groups = x$group,
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


