#' Get the inversion status of a sample
#'
#' This function estimates the inversion status of the samples using the similarity scores and a LDA
#' model.
#'
#' @param scores Matrix of similarity scores (from classifSNPspar)
#' to a group is smaller than this parameter, the individual is classified as NA.
#' @return Factor with the inversion classification of the samples
getInvStatus <- function(scores) {
  class <- factor(colnames(scores)[max.col(scores)], levels = colnames(scores))
  names(class) <- rownames(scores)
  posterior <- scores/rowSums(scores)
  certainty <- apply(posterior, 1, max)
  return(list(class = class, certainty = certainty))
}
