library(adegenet)
library(vcfR)
library(ape)
library(knitr)
library(reshape2)
library(ggplot2)

vcf_file <- "CBS14422_vs_Pc2113_140Mb.F2.recode.vcf"

vcf <- read.vcfR(vcf_file, verbose = FALSE )

### Ploidy I ###

# This code was adopted from the tutorial below: 
# https://knausb.github.io/vcfR_documentation/determining_ploidy_1.html

# The ‘AD’ field in our VCF data includes the depth at which each allele was sequenced.
ad <- extract.gt(vcf, element = 'AD')
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)
ad <- extract.gt(vcf, element = 'AD')
hist(ad2[,"CBS14422"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"CBS14422"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))

#Make these plots look better by looking at allele depth 
ad <- extract.gt(vcf, element = 'AD')
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)

# Subset to a vector for manipulation.
tmp3 <- allele1[,"CBS14422"]
tmp3 <- tmp3[tmp3 <= 100]
hist(tmp3, breaks=seq(0,250,by=1), col="#808080", main = "CBS14422")
sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.1, 0.90), na.rm=TRUE)
sums[,"CBS14422"]
abline(v=sums[,"CBS14422"], col=2, lwd=2)

# We can use this information to filter our data. 
# First we’ll look at the raw data.
boxplot(allele1, las=3)
sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.1, 0.90), na.rm=TRUE)
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

# Now we can check our work with another set of box and whisker plots.
ad <- extract.gt(vcf, element = 'AD')
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
boxplot(allele1, las=3)

# Now we can see if our histogram of allele balance has been cleaned up.
gt <- extract.gt(vcf, element = 'GT')
hets <- is_het(gt)
is.na( ad[ !hets ] ) <- TRUE
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)
# CBS144.22
hist(ad2[,"CBS14422"], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
hist(ad1[,"CBS14422"], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))
