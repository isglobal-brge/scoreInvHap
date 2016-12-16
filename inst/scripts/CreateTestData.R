# Create test data
## Load
library(brgeDummyData)
data("snpsVCF")

testdata <- snpsVCF

## Select SNPs of the VCF that are present in the inversion
invSNPs <- intersect(colnames(snpsVCF$genotypes),  names(SNPsR2$ROIno.8.3))

testdata$genotypes <- testdata$genotypes[, invSNPs[1:]]
