library(GenomicRanges)
a <- read.csv2("~/../Desktop/Libro1.csv", as.is = TRUE)

colnames(a)[c(1, 4:6)] <- c("scoreInvHap.name", "Inv.freq", "Haplotypes", "Num.SNPs")
a[, "Inv.freq"] <- as.numeric(a[, "Inv.freq"])

inversionGR <- GRanges(gsub(",", "", a$Coordinates))
mcols(inversionGR) <- a[, -3]
names(inversionGR) <- a[, 1]
devtools::use_data(inversionGR)
