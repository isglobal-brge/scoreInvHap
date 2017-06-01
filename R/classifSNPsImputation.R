#' @rdname classifSNPs
#' @export
classifSNPsImpute <- function(genos, R2, refs, BPPARAM = BiocParallel::bpparam()){

    # Select SNPs present in R2, references and genotypes
    common <- Reduce(intersect, list(names(R2), names(refs), rownames(genos)))
    R2 <- R2[common]
    genos <- genos[common, , , drop = FALSE]
    refs <- refs[common]
    numrefs <- nrow(refs[[1]])

    # Compute the scores and the probabilities
    scores <- computeScoreImpute(genos, refs = refs, R2 = R2, BPPARAM = BPPARAM)
    snps <- rep(length(common), ncol(genos))
    names(snps) <- names(scores)

    list(scores = scores, numSNPs = snps)
}


computeScoreImpute <- function(geno, refs, R2, BPPARAM = BiocParallel::bpparam()){

    score_list <- BiocParallel::bplapply(rownames(geno), function(snp)
        geno[snp, , ] %*% t(refs[[snp]])*R2[snp], BPPARAM = BPPARAM)
    mat <- Reduce(`+`, score_list)
    mat <- mat/sum(R2[rownames(geno)])
    mat
}
