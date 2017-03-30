#' Get similarity scores
#'
#' This function computes the similarity scores between the sample SNPs and the inversion's reference.
#'
#' @details This function computes two similarity scores for each individual: the inverted and the
#' standard score. The individuals' score are computed as a weighted mean of the similarity scores
#' of each SNP.
#'
#' For each SNP, we compare the genotype of the invidual with the genotypes frequencies
#' in the homozygous inverted and standard samples. If the frequency of the individual genotype is
#' higher in the inverted samples, we add 1 to the inverted score and 0 to the standard score. If the
#' frequency of the individual genotype is higher in the standard sample, we add 1 to the standard
#' score and 0 to the inverted score. The individuals' scores are computed using the mean of the
#' SNPs scores weigthed by the R2 between the SNP and the inversion status.
#'
#' @export
#'
#' @param genos Matrix with the samples genotypes. It is the result of \code{getGenotypesTable}
#' @param R2 Vector with the R2 between the SNPs and the inversion status
#' @param refs List of matrices. Each matrix has, for an SNP, the frequencies of each genotype in the
#' homozygous inverted, homozygous standard and heterozygous samples.
#' @param mc.cores Numeric with the number of cores used in the computation
#' @return Matrix with the standard and inverted scores for each individual. Samples are in columns and
#' scores in rows
#'
classifSNPspar <- function(genos, R2, refs, alfreq, genofreq, mc.cores){

    # Select SNPs present in R2, references and genotypes
    common <- Reduce(intersect, list(names(R2), names(refs), colnames(genos), names(alfreq)))
    R2 <- R2[common]
    genos <- genos[, common, drop = FALSE]
    refs <- refs[common]
    alfreq <- alfreq[common]
    numrefs <- nrow(refs[[1]])

    res <-  parallel::mclapply(rownames(genos), function(ind) {
      computeScore(genos[ind, ], refs = refs, R2 = R2, alfreq = alfreq, genofreq = genofreq)
    }, mc.cores = mc.cores)
    names(res) <- rownames(genos)

    scores <- t(vapply(res, `[[`, numeric(numrefs), "score"))
    probs <- t(vapply(res, `[[`, numeric(numrefs), "prob"))
    snps <- vapply(res, `[[`, numeric(1), "numSNPs")

    list(scores = scores, probs = probs, numSNPs = snps)

    # res <- matrix(unlist(res), nrow = numrefs)
    # rownames(res) <- rownames(refs[[1]])
    # colnames(res) <- rownames(genos)
    # sum_R2 <- matrix(apply(genos != "NN", 1, function(x) sum(R2[x])),
    #                   ncol = nrow(genos), nrow = nrow(res), byrow = T)
    # numSNPs <- rowSums(genos != "NN")
    # scores <- res/sum_R2
    # list(scores = scores, numSNPs = numSNPs)
  }


computeScore <- function(geno, refs, R2, alfreq, genofreq){
  goodgenos <- geno != "NN"
  geno <- geno[goodgenos]
  refs <- refs[goodgenos]
  R2 <- R2[goodgenos]
  numSNPs <- sum(goodgenos)

  mat <- t(sapply(names(geno), function(snp) {
    refs[[snp]][, geno[snp]]*R2[snp]
  }))
  score <- colSums(mat)/sum(R2)

  mat <- log10(t(sapply(names(geno), function(snp) {
    refs[[snp]][, geno[snp]]
  })))
  problog <- colSums(mat)

  haplofreqlog <- sum(log10(sapply(names(geno), function(snp) {
    alfreq[[snp]][geno[snp]]
  })))

  postprob <- 10^(problog-haplofreqlog)*genofreq


  list(score = score, numSNPs = numSNPs, prob = postprob)
}
