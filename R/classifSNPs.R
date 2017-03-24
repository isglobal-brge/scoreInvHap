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
classifSNPspar <- function(genos, R2, refs, mc.cores){

    # Select SNPs present in R2, references and genotypes
    common <- Reduce(intersect, list(names(R2), names(refs), colnames(genos)))
    R2 <- R2[common]
    genos <- genos[, common, drop = FALSE]
    refs <- refs[common]
    numrefs <- nrow(refs[[1]])

    if (length(common) == 0){
      res <- matrix(0, nrow = 2, ncol = nrow(genos))
      rownames(res) <- c("inv", "std")
      colnames(res) <- rownames(genos)
      return(res)
    }

    res <-  parallel::mclapply(rownames(genos),
                     function(ind) {
                       rowSums(sapply(colnames(genos), function(geno){
                         gen <- genos[ind, geno]
                         a <- rep(0, numrefs)
                         if (gen %in% colnames(refs[[geno]])){
                           a <- refs[[geno]][, gen]*R2[geno]
                         }
                         a
                       }))
                     }, mc.cores = mc.cores)
    res <- matrix(unlist(res), nrow = numrefs)
    rownames(res) <- rownames(refs[[1]])
    colnames(res) <- rownames(genos)
    res <- res/sum(R2[colnames(genos)])
    res
  }
