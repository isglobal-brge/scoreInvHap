#' scoreInvHapRes instances
#'
#' Container with the results of the classification pipeline
#'
#' @export
#' @rdname scoreInvHapRes-class
#' @name scoreInvHapRes
#' @aliases scoreInvHapRes-class scoreInvHapRes-methods
#'
#' @slot classification Factor with the individuals classification
#' @slot scores Simmilarity scores for the different haplotypes.
#' @slot numSNPs Numeric with SNPs used to compute the scores.
#' @slot certainty Numeric with the certainty of the classification for each individual.
#' @return A scoreInvHapRes instance
#' @examples
#' if(require(VariantAnnotation)){
#'     vcf <- readVcf(system.file("extdata", "example.vcf", package = "scoreInvHap"), "hg19")
#'
#'     ## Create scoreInvHapRes class from pipeline
#'     res <- scoreInvHap(vcf, inv = inv7_005)
#'
#'     ## Print object
#'     res
#'
#'     ## Get haplotype classification
#'     classification(res)
#'
#'     ## Get similiraty scores
#'     scores(res)
#' }
setClass (
    Class = "scoreInvHapRes",
    representation(
        classification = "factor",
        scores = "matrix",
        numSNPs = "numeric",
        certainty = "numeric"
    )
)

#' @export
setGeneric("classification", function(object, minDiff = 0, callRate = 0, inversion = FALSE){
    standardGeneric("classification")
})

#' @describeIn scoreInvHapRes Get classification
#' @aliases scoreInvHapRes-methods classification
#' @param object \code{scoreInvHapRes}
#' @param minDiff Numeric with the threshold of the minimum difference
#' between the top and the second score. Used to filter samples.
#' @param callRate Numeric with the threshold of the minimum call rate
#' of the samples. Used to filter samples.
#' @param inversion Logical. If true, haplotypes classification is adapted
#' to return inversion status.
setMethod(
    f = "classification",
    signature = "scoreInvHapRes",
    definition = function(object, minDiff, callRate, inversion) {
        res <- object@classification
        goodScores <- diffscores(object) > minDiff
        goodSNPs <- propSNPs(object) > callRate

        res <- res[goodScores & goodSNPs]
        if (inversion){
            regmatches(levels(res), gregexpr("Na", levels(res)), invert = FALSE) <- "N"
            regmatches(levels(res), gregexpr("Nb", levels(res)), invert = FALSE) <- "N"
            regmatches(levels(res), gregexpr("Nc", levels(res)), invert = FALSE) <- "N"
            regmatches(levels(res), gregexpr("Ia", levels(res)), invert = FALSE) <- "I"
            regmatches(levels(res), gregexpr("Ib", levels(res)), invert = FALSE) <- "I"
            regmatches(levels(res), gregexpr("Ic", levels(res)), invert = FALSE) <- "I"
        }
        return(res)
    }
)


#' @export
setGeneric("certainty", function(object){
    standardGeneric("certainty")
})

#' @describeIn scoreInvHapRes Get classification certainty
#' @aliases scoreInvHapRes-methods certainty
setMethod(
    f = "certainty",
    signature = "scoreInvHapRes",
    definition = function(object) {
        return(object@certainty)
    }
)

#' @export
setGeneric("diffscores", function(object){
    standardGeneric("diffscores")
})

#' @describeIn scoreInvHapRes Get maximum similarity scores
#' @aliases scoreInvHapRes-methods diffscores
setMethod(
    f = "diffscores",
    signature = "scoreInvHapRes",
    definition = function(object) {
        return(apply(scores(object), 1, function(x) max(x) - sort(x, decreasing = TRUE)[2]))
    }
)

#' @export
setGeneric("maxscores", function(object){
    standardGeneric("maxscores")
})

#' @describeIn scoreInvHapRes Get maximum similarity scores
#' @aliases scoreInvHapRes-methods maxscores
setMethod(
    f = "maxscores",
    signature = "scoreInvHapRes",
    definition = function(object) {
        return(apply(scores(object), 1, max))
    }
)

#' @export
setGeneric("numSNPs", function(object){
    standardGeneric("numSNPs")
})

#' @describeIn scoreInvHapRes Get number of SNPs used in computation
#' @aliases scoreInvHapRes-methods numSNPs
setMethod(
    f = "numSNPs",
    signature = "scoreInvHapRes",
    definition = function(object) {
        return(object@numSNPs)
    }
)

#' @export
setGeneric("plotCallRate", function(object, callRate = 0.9, ...){
    standardGeneric("plotCallRate")
})

#' @describeIn scoreInvHapRes Plot call rate based QC
#' @aliases scoreInvHapRes-methods plotCallRate
#' @param ... Further parameters passed to plot function.
setMethod(
    f = "plotCallRate",
    signature = "scoreInvHapRes",
    definition = function(object, callRate = 0.9, ...) {

        graphics::hist(propSNPs(object), breaks = seq(0, 1, 0.05), xlim = c(0, 1),
             xlab = "SNPs proportion", ...)
        graphics::abline(v = callRate)
    }
)

#' @export
setGeneric("plotScores", function(object, minDiff = 0.1, ...){
    standardGeneric("plotScores")
})

#' @describeIn scoreInvHapRes Plot scores based QC
#' @aliases scoreInvHapRes-methods plotScores
setMethod(
    f = "plotScores",
    signature = "scoreInvHapRes",
    definition = function(object, minDiff = 0.1, ...) {

        graphics::plot(maxscores(object), diffscores(object), xlim = c(0, 1),
             xlab = "Max Scores", ylab = "Diff Score", ...)
        graphics::abline(h = minDiff)
    }
)


#' @export
setGeneric("propSNPs", function(object){
    standardGeneric("propSNPs")
})

#' @export
#' @describeIn scoreInvHapRes Get proportions of SNPs used in computation
#' @aliases scoreInvHapRes-methods propSNPs
setMethod(
    f = "propSNPs",
    signature = "scoreInvHapRes",
    definition = function(object) {
        return(numSNPs(object)/max(numSNPs(object)))
    }
)

#' @export
setGeneric("scores", function(object){
    standardGeneric("scores")
})

#' @describeIn scoreInvHapRes Get similarity scores
#' @aliases scoreInvHapRes-methods scores
setMethod(
    f = "scores",
    signature = "scoreInvHapRes",
    definition = function(object) {
        return(object@scores)
    }
)

setMethod(
    f = "show",
    signature = "scoreInvHapRes",
    definition = function(object) {
        class <- classification(object)
        cat("scoreInvHapRes\n")
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
        # cat(sprintf("Mean certainty: %.4f\n", mean(certainty(object), na.rm = TRUE)))

    }
)
