#' Compute the allele table
#'
#' Get a data.frame that maps the numeric genotype of a SNPmatrix (0, 1, 2) into the real genotype.
#' Heterozygous genotypes are ordered alphabetically.
#'
#' @param map Data.frame with the annotation of the SNPs (from plink format)
#' @return Data.frame with genotypes map
getAlleleTable <- function(map){
    AB <- matrix(c(as.character(map$allele.1), as.character(map$allele.2)), ncol = 2)
    allele <- data.frame(AA = paste0(map$allele.1, map$allele.1),
                         AB = apply(AB, 1, function(x) paste(sort(x), collapse = "")),
                         BB = paste0(map$allele.2, map$allele.2), stringsAsFactors = FALSE)
    rownames(allele) <- rownames(map)
    allele$Fail <- rep("NN", nrow(allele))
    allele
}
