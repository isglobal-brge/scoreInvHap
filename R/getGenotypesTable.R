#' Get genotypes table
#'
#' Get a matrix with the sample genotypes from all SNP.
#'
#' @param geno SnpMatrix (from plink format)
#' @param allele Data.frame with the alleles per SNP (from getAlleleTable)
#' @return Character matrix with the samples genotypes
getGenotypesTable <- function(geno, allele){
  names <- rownames(geno)

  ## Convert matrix of genotypes to numeric (0, 1, 2). Add 1 to match column indeces of allele table.
  geno <- methods::as(geno, "numeric") + 1
  geno[is.na(geno)] <- 4
  geno <- vapply(1:ncol(geno), function(x) unlist(allele[x, geno[,x], drop = TRUE]),
                 character(nrow(geno)))
  colnames(geno) <- rownames(allele)
  rownames(geno) <- names
  geno
}
