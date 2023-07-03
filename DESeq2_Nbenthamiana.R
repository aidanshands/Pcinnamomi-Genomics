# Aidan Shands 

# R v.4.0.3 
# DEseq2 v.1.30.1
# import libraries (may not need all)
library(DESeq2)
library(LSD)
library(genefilter)
library(biomaRt)
library(dplyr)
library(limma)
#-------------------------------------------------------------------------------
# I followed this guide in 2021:
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# by Michael I. Love, Simon Anders, and Wolfgang Huber
#-------------------------------------------------------------------------------
# Importing Count data 
countData <- as.matrix(read.csv("GX_Nben_PcSomMyc_Count_Matix.csv", header = T, 
                                fill = T, row.names = 1))

#Setting up Time for ColData 
# Time 0 represents the 3 P. cinnamomi somatic tissue samples
Time = c(rep("0h",3),rep("16h",3),rep("24h",3))
colData <- data.frame(row.names = colnames(countData), Time = Time)

# Setting up DEseq count matrix 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, 
                              ignoreRank = FALSE, design = ~Time)

# minimal pre-filtering to keep only rows that have at least 10 reads total.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#-------------------------------------------------------------------------------
# Differential Expression Analysis
#-------------------------------------------------------------------------------
# The "poscounts" estimator deals with a gene with some zeros, 
# by calculating a modified geometric mean by taking the n-th root of the 
# product of the non-zero counts. 
# See https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf
dds <- estimateSizeFactors(dds, type = 'poscounts')
dds <- DESeq(dds, test="Wald")
dds_results <- results(dds)

# Normalized counts
N_Counts = counts(dds, normalized=T)
# write.csv(N_Counts, file = "Nb_Normalized_counts.csv")

#-------------------------------------------------------------------------------
# Comparisons
#-------------------------------------------------------------------------------
# 16h vs Mycelium
res_0_vs_16 <- results(dds, contrast=c("Time","16h","0h"), 
                       alpha = 0.05, lfcThreshold=1)
# write.csv(res_0_vs_16, file = "Nben_16h.csv")
summary(res_0_vs_16)

# 24h vs Mycelium
res_0_vs_24 <- results(dds, contrast=c("Time","24h","0h"), 
                          alpha = 0.05, lfcThreshold=1)
# write.csv(res_0_vs_24, file = "Nben_24h.csv")
summary(res_0_vs_24)