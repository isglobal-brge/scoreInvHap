context("Test scoreInvHap")
library(VariantAnnotation)
vcf_file <- system.file("extdata", "example.vcf", package = "scoreInvHap")
vcf <- readVcf(vcf_file, "hg19")

res1 <- scoreInvHap(SNPlist = vcf, inv = "inv7_005")
res1

res2 <- scoreInvHap(SNPlist = vcf, inv = "inv7_005", probs=TRUE)
res2

# microbenchmark::microbenchmark(res1 = scoreInvHap(SNPlist = vcf,
#                                                  inv = "inv7_005"),
#                               res2 = scoreInvHap(SNPlist = vcf,
#                                                  imputed=TRUE,
#                                                  inv = "inv7_005"))

