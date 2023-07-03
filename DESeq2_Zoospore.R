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
countData <- as.matrix(read.csv("Zoospore_Count_Matrix.csv", 
                                header = T, fill = T, row.names = 1))


#Setting up Time for ColData 
Condition = c(rep("Myc",3),rep("Zoo",3))
colData <- data.frame(row.names = colnames(countData), Condition = Condition)


# Setting up DEseq count matrix 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, 
                              ignoreRank = FALSE, design = ~Condition)

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
#write.csv(N_Counts, file = "Zoo_Normalized_counts.csv")

#-------------------------------------------------------------------------------
# Comparisons
#-------------------------------------------------------------------------------

Zoo_vs_Myc <- results(dds, contrast=c("Condition","Zoo","Myc"), 
                      alpha = 0.05, lfcThreshold=1)
# write.csv(Zoo_vs_Myc, file = "Zoospore.csv")
summary(Zoo_vs_Myc_v2)
