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
#' @slot certainty Numeric with the certainty of the classification for each individual.
#' @slot scores Simmilarity scores for the inverted and standard references.
setClass (
  Class = "SNPfieRes",
  representation(
    classification = "factor",
    certainty = "numeric",
    scores = "matrix"
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
#' @param object \code{SNPfieRes}
setMethod(
  f = "certainty",
  signature = "SNPfieRes",
  definition = function(object) {
    return(object@certainty)
  }
)


#' @export
setGeneric("scores", function(object){
  standardGeneric("scores")
})

#' @describeIn SNPfieRes Get samples similarity scores
#' @aliases SNPfieRes-methods scores
#' @param object \code{SNPfieRes}
setMethod(
  f = "scores",
  signature = "SNPfieRes",
  definition = function(object) {
    return(object@scores)
  }
)

#' @describeIn SNPfieRes Plot samples similarity scores
#' @aliases SNPfieRes-methods plot
#' @param object \code{SNPfieRes}
setMethod(
  f = "plot",
  signature = "SNPfieRes",
  definition = function(x, ...) {
    sc <- scores(x)
    cl <- classification(x)
    colors <- c("red", "green", "blue")
    names(cl) <- c("NI/NI", "NI/I", ("I/I"))
    plot(sc[, 1], sc[, 2], col = colors[cl], xlab = "Standard Score",
         ylab = "Inverted score", xlim = c(0, 1), ylim = c(0, 1), ...)
    legend("topright", c("NI/NI", "NI/I", "I/I"), pch = 16, col = colors)
  }
)

setMethod(
  f = "show",
  signature = "SNPfieRes",
  definition = function(object) {
    class <- classification(object)
    cat("SNPfieRes\n")
    cat("Samples: ", length(class), "\n")
    cat("Genotypes' table:\nNN\tNI\tII\n", sum(class == "NI/NI"), "\t",
        sum(class == "NI/I"), "\t", sum(class == "I/I"), "\n")
    cat(sprintf("Inversion frequency: %.2f%%\n",
            sum(c(class == "NI/I", (class == "I/I") * 2))*100/(length(class)*2)))
    cat(sprintf("Mean certainty: %.4f\n", mean(certainty(object))))

  }
)
