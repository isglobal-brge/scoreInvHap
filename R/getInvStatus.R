#' Get the inversion status of a sample
#'
#' This function estimates the inversion status of the samples using the similarity scores and a LDA
#' model.
#'
#' @param scores Matrix of similarity scores (from classifSNPspar)
#' @param ldamodel LDAmodel used to classify individuals
#' @param uncertainy.thres Numeric with the uncertainty threshold. If the maximum certainty of belonging
#' to a group is smaller than this parameter, the individual is classified as NA.
#' @return Factor with the inversion classification of the samples
getInvStatus <- function(scores, ldamodel) {
  pred <- MASS:::predict.lda(ldamodel, data.frame(t(scores)))
  class <- pred$class
  names(class) <- colnames(scores)
  certainty <- apply(pred$posterior, 1, max)
  res <- list(class = class, certainty = certainty)
  return(res)
}
