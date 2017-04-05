#' Get the inversion status of a sample
#'
#' This function estimates the inversion status of the samples using the probabilities
#' computed in \code{classifSNPs}
#'
#' @param scores Matrix of probabilities (from \code{classifSNPs})
#' @return List with the results:
#' \itemize{
#' \item{class: Vector with the most probable classification}
#' \item{certainty: Vector with the certainty of the most probable classification}
#' }
getInvStatus <- function(scores) {
  scores[scores == Inf] <- 1
  class <- factor(colnames(scores)[max.col(scores)], levels = colnames(scores))
  names(class) <- rownames(scores)
  posterior <- scores/rowSums(scores)
  certainty <- apply(posterior, 1, max)
  return(list(class = class, certainty = certainty))
}
