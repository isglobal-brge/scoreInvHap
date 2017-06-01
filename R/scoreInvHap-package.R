#' scoreInvHap: package to get inversion status of predefined regions.
#'
#' scoreInvHap can get the samples' inversion status of known inversions. scoreInvHap uses SNP data as input
#' and requires the following information about the inversion: genotype frequencies in the different
#' inversion groups, R2 between the region SNPs and inversion status, heterozygote genotypes in the
#' reference, allele frequencies in the reference population and inversion frequencies.
#' The package include this data for two well known inversions (8p23 and
#' 17q21.31) and for two additional validated regions.
#'
#' @docType package
#' @name scoreInvHap
#'
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom Biostrings complement DNAStringSet
#' @importFrom GenomicRanges mcols
#' @importFrom methods as
#' @importClassesFrom snpStats SnpMatrix
NULL
