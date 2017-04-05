#' SNPfieRes instances
#'
#' Container with the results of the classification pipeline
#'
#' @export
#' @rdname SNPfieRes-class
#' @name SNPfieRes
#' @aliases SNPfieRes-class SNPfieRes-methods
#'
#' @slot classification Factor with the individuals classification
#' @slot probs Probability of each individual to belong to the different haplotypes
#' @slot scores Simmilarity scores for the different haplotypes.
#' @slot numSNPs Numeric with SNPs used to compute the scores.
#' @slot certainty Numeric with the certainty of the classification for each individual.
setClass (
  Class = "SNPfieRes",
  representation(
    classification = "factor",
    probs = "matrix",
    scores = "matrix",
    numSNPs = "numeric",
    certainty = "numeric"
  )
)

#' @export
setGeneric("classification", function(object){
  standardGeneric("classification")
})

#' @describeIn SNPfieRes Get samples classification
#' @aliases SNPfieRes-methods classification
#' @param object \code{SNPfieRes}
setMethod(
  f = "classification",
  signature = "SNPfieRes",
  definition = function(object) {
    return(object@classification)
  }
)


#' @export
setGeneric("certainty", function(object){
  standardGeneric("certainty")
})

#' @describeIn SNPfieRes Get classification certainty
#' @aliases SNPfieRes-methods certainty
setMethod(
  f = "certainty",
  signature = "SNPfieRes",
  definition = function(object) {
    return(object@certainty)
  }
)

#' @export
setGeneric("numSNPs", function(object){
  standardGeneric("numSNPs")
})

#' @describeIn SNPfieRes Get number of SNPs used in computation
#' @aliases SNPfieRes-methods numSNPs
setMethod(
  f = "numSNPs",
  signature = "SNPfieRes",
  definition = function(object) {
    return(object@numSNPs)
  }
)

#' @export
setGeneric("probs", function(object){
  standardGeneric("probs")
})

#' @describeIn SNPfieRes Get samples probabilities
#' @aliases SNPfieRes-methods probs
setMethod(
  f = "probs",
  signature = "SNPfieRes",
  definition = function(object) {
    return(object@probs)
  }
)


#' @export
setGeneric("scores", function(object){
  standardGeneric("scores")
})

#' @describeIn SNPfieRes Get samples similarity scores
#' @aliases SNPfieRes-methods scores
setMethod(
  f = "scores",
  signature = "SNPfieRes",
  definition = function(object) {
    return(object@scores)
  }
)

# #' @describeIn SNPfieRes Plot samples similarity scores
# #' @aliases SNPfieRes-methods plot
# #' @param object \code{SNPfieRes}
# setMethod(
#   f = "plot",
#   signature = "SNPfieRes",
#   definition = function(x, ...) {
#     sc <- scores(x)
#     cl <- as.character(classification(x))
#     colors <- c("red", "green", "blue")
#     names(colors) <- c("NI/NI", "NI/I", ("I/I"))
#     plot(sc[, 1], sc[, 2], col = colors[cl], xlab = "Standard Score",
#          ylab = "Inverted score", xlim = c(0, 1), ylim = c(0, 1), ...)
#     legend("topright", c("NI/NI", "NI/I", "I/I"), pch = 16, col = colors)
#   }
# )

setMethod(
  f = "show",
  signature = "SNPfieRes",
  definition = function(object) {
    class <- classification(object)
    cat("SNPfieRes\n")
    cat("Samples: ", length(class), "\n")
    tab <- table(class)
    cat("Genotypes' table:\n", paste0(names(tab), "\t"), "\n", paste0(tab, "\t"), "\n")
    if (all(grepl("N", names(tab)) | grepl("I", names(tab)))){
      N <- lengths(regmatches(names(tab), gregexpr("N", names(tab))))
      invtab <- c(NN = sum(tab[N == 2]), NI = sum(tab[N == 1]),
                  II = sum(tab[N == 0]))
    cat("- Inversion genotypes' table:\n",
        paste0(names(invtab), "\t"), "\n", paste0(invtab, "\t"), "\n")
    cat(sprintf("- Inversion frequency: %.2f%%\n",
                sum(invtab[2], 2*invtab[3])/(sum(invtab)*2)*100))
    }
    cat(sprintf("Mean certainty: %.4f\n", mean(certainty(object), na.rm = TRUE)))

  }
)
