#' SNPfier: package to get inversion status of predefined regions.
#'
#' SNPfier can get the samples' inversion status of known inversions. SNPfier uses SNP data as input
#' and requires the following information about the inversion: genotype frequencies in the different
#' inversion groups, R2 between the region SNPs and inversion status, heterozygote genotypes in the
#' reference, allele frequencies in the reference population and inversion frequencies.
#' The package include this data for two well known inversions (8p23 and
#' 17q21.31) and for two additional validated regions.
#'
#' @docType package
#' @name SNPfier
#'
#' @import MASS
#' @importFrom Biostrings complement DNAStringSet
#' @importFrom methods as
#' @importFrom parallel mclapply
NULL
