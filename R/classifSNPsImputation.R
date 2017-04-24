#' Get similarity scores and probability using imputed data.
#'
#' This function computes the similarity scores between the sample SNPs and the haplotype's
#' reference. It also returns the probability of each sample to belong to the
#' different haplotypes.
#'
#' @details This function computes, for each individual, similarity scores for
#' all the present haplotypes. For each SNP, we compute as many similarity scores
#' as haplotypes present in the reference. We have defined the similarity score as
#' the frequency of this genotype in the different haplotype population. To compute
#' the global similarity score, we have computed a mean of the scores by SNP weighted
#' by the R2 between the SNP and the haplotype classification.
#'
#' In addition, we have computed the probability of a sample to belong to each of
#' the haplotypes. To do so, we compute the conditional probability of a sample genotype
#' to belong to each haplotype and we apply Bayes' theorem.
#'
#' @export
#'
#' @param genos Array with the genotypes probabilities.
#' @param R2 Vector with the R2 between the SNPs and the inversion status
#' @param refs List of matrices. Each matrix has, for an SNP, the frequencies of each genotype in the
#' different haplotypes.
#' @param mc.cores Numeric with the number of cores used in the computation
#' @return List with the results:
#' \itemize{
#' \item{scores: Matrix with the simmilarity scores of the individuals}
#' \item{numSNPs: Vector with the number of SNPs used in each computation}
#' }
#'
classifSNPsImpute <- function(genos, R2, refs, mc.cores = 1){

  # Select SNPs present in R2, references and genotypes
  common <- Reduce(intersect, list(names(R2), names(refs), rownames(genos)))
  R2 <- R2[common]
  genos <- genos[common, , , drop = FALSE]
  refs <- refs[common]
  numrefs <- nrow(refs[[1]])

  # Compute the scores and the probabilities
  scores <- computeScoreImpute(genos, refs = refs, R2 = R2, mc.cores = mc.cores)
  snps <- rep(length(common), ncol(genos))
  names(snps) <- names(scores)

  list(scores = scores, numSNPs = snps)
}


computeScoreImpute <- function(geno, refs, R2, mc.cores = 1){

  score_list <- mclapply(rownames(geno), function(snp)
    geno[snp, , ] %*% t(refs[[snp]])*R2[snp], mc.cores = mc.cores)
  mat <- Reduce(`+`, score_list)
  mat <- mat/sum(R2[rownames(geno)])
  mat
}
