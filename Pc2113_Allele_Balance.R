library(adegenet)
library(vcfR)
library(ape)
library(knitr)
library(reshape2)
library(ggplot2)

vcf_file <- "Pc2113_Pc2109_vs_Pc2113T1.F2.recode.vcf"

vcf <- read.vcfR(vcf_file, verbose = FALSE )

### Ploidy I ###
# https://knausb.github.io/vcfR_documentation/determining_ploidy_1.html
ad <- extract.gt(vcf, element = 'AD')
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)
ad <- extract.gt(vcf, element = 'AD')
hist(ad2[,"Pc2113"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"Pc2113"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
ad <- extract.gt(vcf, element = 'AD')
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
tmp3 <- allele1[,"Pc2113"]
tmp3 <- tmp3[tmp3 <= 250]
hist(tmp3, breaks=seq(0,250,by=1), col="#808080", main = "Pc2113")
sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.075, 0.95), na.rm=TRUE)
sums[,"Pc2113"]
abline(v=sums[,"Pc2113"], col=2, lwd=2)
boxplot(allele1, las=3)
sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.075, 0.95), na.rm=TRUE)
# Allele 1
dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[1,])
vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[2,])
vcf@gt[,-1][dp2 > 0] <- NA
# Allele 2
dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[1,])
vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[2,])
vcf@gt[,-1][dp2 > 0] <- NA
ad <- extract.gt(vcf, element = 'AD')
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
boxplot(allele1, las=3)
gt <- extract.gt(vcf, element = 'GT')
hets <- is_het(gt)
is.na( ad[ !hets ] ) <- TRUE
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)
hist(ad2[,"Pc2113"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"Pc2113"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))