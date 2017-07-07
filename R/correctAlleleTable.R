#' Solve genotypes discrepancies
#'
#' This function tries to solve discrepancies between the reference and sample genotypes. The cause of
#' these discrepancies is that samples and references have used different strands to codify the SNP. This
#' function get the complement genotypes for the discordant SNPs and checks if discordancies are solved.
#'
#' @param alleletable Data.frame with the alleles per SNP (from getAlleleTable)
#' @param hetRefs Character vector with the heterozygous genotypes in the reference.
#' @param map Data.frame with the annotation of the SNPs (from plink format)
#' @return alleletable without discrepancies between these genotypes and the references.
correctAlleleTable <- function(alleletable, hetRefs, map){
    alleletable$InvRef <- unlist(hetRefs)[rownames(alleletable)]

    # Select SNPs that the alleles on the array are not in our reference
    wrongSNPs <- rownames(alleletable)[sapply(rownames(alleletable),
                                              function(x) !alleletable[x, 5] %in% alleletable[x, 1:3])]

    badMap <- map[wrongSNPs, ]
    badMap$allele.1 <- as.character(Biostrings::complement(Biostrings::DNAStringSet(badMap[, "allele.1"])))
    badMap$allele.2 <- as.character(Biostrings::complement(Biostrings::DNAStringSet(badMap[, "allele.2"])))
    alleletable[wrongSNPs, ] <- getAlleleTable(badMap)
    alleletable$InvRef <- unlist(hetRefs)[rownames(alleletable)]
    wrongSNPs <- rownames(alleletable)[sapply(rownames(alleletable),
                                              function(x) !alleletable[x, 5] %in% alleletable[x, 1:3])]
    alleletable <- alleletable[!rownames(alleletable) %in% wrongSNPs, ]
    alleletable
}
