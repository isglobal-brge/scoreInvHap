#' @rdname classifSNPs
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
