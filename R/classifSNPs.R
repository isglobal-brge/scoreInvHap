#' Get similarity scores and probability
#'
#' This function computes the similarity scores between the sample SNPs and the haplotype's
#' reference.
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
#' @param R2 Vector with the R2 between the SNPs and the inversion status.
#' @param refs List of matrices. Each matrix has, for an SNP, the frequencies of each genotype in the
#' different haplotypes.
#' @param alleletable Data frame with the reference alleles computed with \code{getAlleleTable}.
#' @param BPPARAM A \code{BiocParallelParam} instance. Used to parallelize computation
#' @return List with the results:
#' \itemize{
#' \item{scores: Matrix with the simmilarity scores of the individuals}
#' \item{numSNPs: Vector with the number of SNPs used in each computation}
#' }
#' @examples
#' ## Simulate a table of genotypes from inv8_001
#' geno <- matrix(c("CC", "GG", "AA", "CG", "NN", "AC", "GG", "AA", "CC"),
#'     nrow = 3, dimnames = list(letters[1:3],
#'     c("rs141039449", "rs138092889", "rs138217047")))
#'
#' ## Run function using reference of inv8p23.1
#' classifSNPs(geno, SNPsR2$inv8_001, Refs$inv8_001)
classifSNPs <- function(genos, R2, refs, alleletable,
                        BPPARAM = BiocParallel::SerialParam()){

    # Select SNPs present in R2, references and genotypes
    common <- Reduce(intersect, list(names(R2), names(refs), colnames(genos)))
    R2 <- R2[common]
    genos <- genos[, common, drop = FALSE]
    refs <- refs[common]
    numrefs <- nrow(refs[[1]])

    # Compute the scores

    ff <- function(ind, alleletable) {
        genos <- getGenotypesTable(geno = genos[ind, , drop=FALSE],
                                   allele = alleletable)
        computeScore(genos, refs = refs, R2 = R2)
    }

    res <-  BiocParallel::bplapply(rownames(genos), ff,
                                   alleletable=alleletable,
                                   BPPARAM = BPPARAM)
    names(res) <- rownames(genos)

    scores <- t(vapply(res, `[[`, numeric(numrefs), "score"))
    snps <- vapply(res, `[[`, numeric(1), "numSNPs")

    list(scores = scores, numSNPs = snps)

}

#' Compute all similarity scores for a sample
#'
#' @description Internal
#'
#' @param geno Vector with the sample genotypes. It is the result of
#' \code{getGenotypesTable}
#' @param refs List of matrices. Each matrix has, for an SNP, the frequencies of each genotype in the
#' different haplotypes.
#' @param R2 Vector with the R2 between the SNPs and the inversion status
#' @return List with the results:
#' \itemize{
#' \item{scores: Vector with the simmilarity scores of the sample}
#' \item{numSNPs: Numeric with the number of SNPs used in the computation}
#' }
computeScore <- function(geno, refs, R2){

    ## Check if geno has names
    if (is.null(names(geno))){
        names(geno) <- names(refs)
    }

    ## Filter SNPs without a calling
    goodgenos <- geno != "NN"

    numSNPs <- sum(goodgenos)
    haplos <- rownames(refs[[1]])
    numhaplos <- length(haplos)

    ## If any SNP has a calling, return 0 for all scores
    if(numSNPs == 0){
        score <- postprob <- rep(0, numhaplos)
        names(score) <- names(postprob) <- haplos
        return(list(score = score, numSNPs = numSNPs, prob = postprob))
    }

    geno <- geno[goodgenos]
    refs <- refs[goodgenos]
    R2 <- R2[goodgenos]

    ## Modify geno and refs to allow vectorized computation
    geno <- paste0(names(geno), geno)
    refs <- t(Reduce(cbind, lapply(names(refs), function(x) {
        colnames(refs[[x]]) <- paste0(x, colnames(refs[[x]]))
        refs[[x]]
        })))

    ## Get the similarity scores
    mat <- data.frame(refs)[geno, ]*matrix(R2, ncol = numhaplos,
                                  nrow = numSNPs)
    mat[is.na(mat)] <- 0
    score <- colSums(mat)/sum(R2)

    list(score = score, numSNPs = numSNPs)
}
