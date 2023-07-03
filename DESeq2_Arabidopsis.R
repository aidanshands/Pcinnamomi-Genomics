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
countData <- as.matrix(read.csv("At_PcSomMyc_Count_Matrix.csv", header = T, 
                                fill = T, row.names = 1))

# Setting up Time for ColData 
# Time 0 represents the 3 P. cinnamomi somatic tissue samples
Time = c(rep("0h",3),rep("16h",3),rep("24h",3), rep("90h",2))
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
# Run DESeq function
dds <- DESeq(dds, test="Wald")
dds_results <- results(dds)

# Normalized counts
N_Counts = counts(dds, normalized=T)

# write.csv(N_Counts, file = "At_Normalized_counts.csv")

#-------------------------------------------------------------------------------
# Comparisons
#-------------------------------------------------------------------------------
# 16h vs Mycelium
res_0_vs_16 <- results(dds, contrast=c("Time","16h","0h"), 
                       alpha = 0.05, lfcThreshold=1)
write.csv(res_0_vs_16, file = "Arabidopsis_16h.csv")
summary(res_0_vs_16)

# 24h vs Mycelium
res_0_vs_24 <- results(dds, contrast=c("Time","24h","0h"), 
                       alpha = 0.05, lfcThreshold=1)
write.csv(res_0_vs_24, file = "Arabidopsis_24h.csv")
summary(res_0_vs_24)

# 90h vs Mycelium
res_0_vs_90 <- results(dds, contrast=c("Time","90h","0h"), 
                       alpha = 0.05, lfcThreshold=1)
write.csv(res_0_vs_90, file = "Arabidopsis_90h.csv")
summary(res_0_vs_90)
