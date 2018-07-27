#' scoreInvHap: package to get inversion status of predefined regions.
#'
#' scoreInvHap can get the samples' inversion status of known inversions. scoreInvHap uses SNP data as input
#' and requires the following information about the inversion: genotype frequencies in the different
#' inversion groups, R2 between the region SNPs and inversion status, heterozygote genotypes in the
#' reference, allele frequencies in the reference population and inversion frequencies.
#' The package include this data for 20 inversions.
#'
#' @docType package
#' @name scoreInvHap
#'
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom Biostrings complement DNAStringSet
#' @importFrom GenomicRanges mcols start
#' @importFrom graphics abline hist plot
#' @importFrom methods as is new
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom VariantAnnotation info
#' @importClassesFrom snpStats SnpMatrix
NULL
