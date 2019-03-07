#' @rdname classifSNPs
#' @export
classifSNPsImpute <- function(genos, R2, refs, BPPARAM = BiocParallel::SerialParam()){

    # Select SNPs present in R2, references and genotypes
    common <- Reduce(intersect, list(names(R2), names(refs), rownames(genos)))
    R2 <- R2[common]
    genos <- genos[common, , , drop = FALSE]
    refs <- refs[common]
    numrefs <- nrow(refs[[1]])

    # Compute the scores and the probabilities
    score_list <- BiocParallel::bplapply(rownames(genos), function(snp)
        tcrossprod(genos[snp, , ],refs[[snp]])*R2[snp],
        BPPARAM = BPPARAM)
    scores <- Reduce(`+`, score_list)
    scores <- scores/sum(R2[rownames(genos)])
    snps <- rep(length(common), ncol(genos))
    names(snps) <- names(scores)

    list(scores = scores, numSNPs = snps)
}
