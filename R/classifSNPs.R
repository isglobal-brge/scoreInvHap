#' Get similarity scores and probability
#'
#' This function computes the similarity scores between the sample SNPs and the haplotype's
#' reference. It also returns the probability of each sample to belong to the
#' different haplotypes.
#'
#' @details classifSNPs computes, for each individual, similarity scores for
#' all the present haplotypes. For each SNP, we compute as many similarity scores
#' as haplotypes present in the reference. We have defined the similarity score as
#' the frequency of this genotype in the different haplotype population. To compute
#' the global similarity score, we have computed a mean of the scores by SNP weighted
#' by the R2 between the SNP and the haplotype classification.
#'
#' classifSNPsImpute is a version of classifSNPs that works with posterior probabilities
#' of imputed genotypes.
#'
#' @export
#'
#' @param genos Matrix with the samples genotypes. It is the result of \code{getGenotypesTable}
#' @param R2 Vector with the R2 between the SNPs and the inversion status
#' @param refs List of matrices. Each matrix has, for an SNP, the frequencies of each genotype in the
#' different haplotypes.
#' @param mc.cores Numeric with the number of cores used in the computation
#' @return List with the results:
#' \itemize{
#' \item{scores: Matrix with the simmilarity scores of the individuals}
#' \item{numSNPs: Vector with the number of SNPs used in each computation}
#' }
classifSNPs <- function(genos, R2, refs, mc.cores){

    # Select SNPs present in R2, references and genotypes
    common <- Reduce(intersect, list(names(R2), names(refs), colnames(genos)))
    R2 <- R2[common]
    genos <- genos[, common, drop = FALSE]
    refs <- refs[common]
    numrefs <- nrow(refs[[1]])

    # Compute the scores and the probabilities
    res <-  parallel::mclapply(rownames(genos), function(ind) {
      computeScore(genos[ind, ], refs = refs, R2 = R2)
    }, mc.cores = mc.cores)
    names(res) <- rownames(genos)

    scores <- t(vapply(res, `[[`, numeric(numrefs), "score"))
    snps <- vapply(res, `[[`, numeric(1), "numSNPs")

    list(scores = scores, numSNPs = snps)

  }


computeScore <- function(geno, refs, R2){

  if (is.null(names(geno))){
    names(geno) <- names(refs)
  }

  goodgenos <- geno != "NN"

  numSNPs <- sum(goodgenos)
  haplos <- rownames(refs[[1]])
  numhaplos <- length(haplos)

  if(numSNPs == 0){
    score <- postprob <- rep(0, numhaplos)
    names(score) <- names(postprob) <- haplos
    return(list(score = score, numSNPs = numSNPs, prob = postprob))
  }

  geno <- geno[goodgenos]
  refs <- refs[goodgenos]
  R2 <- R2[goodgenos]
  mat <- t(vapply(names(geno), function(snp) {
    tryCatch(refs[[snp]][, geno[snp]]*R2[snp], error = function(e) {
      res <- rep(0, numhaplos)
      names(res) <- haplos
      res
      })
  }, numeric(numhaplos)))

  score <- colSums(mat)/sum(R2)


  list(score = score, numSNPs = numSNPs)
}
