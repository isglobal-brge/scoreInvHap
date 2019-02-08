context("Test scoreInvHap")
library(VariantAnnotation)
vcf_file <- system.file("extdata", "example.vcf", package = "scoreInvHap")
vcf <- readVcf(vcf_file, "hg19")

res <- scoreInvHap(SNPlist = vcf, inv = "inv7_005")
res
