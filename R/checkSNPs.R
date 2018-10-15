#' Check genotype object
#'
#' This function checks the genotype object before passing the SNPs to `scoreInvHap`. The
#' function removes SNPs with different alleles or different allele frequencies. Nonetheless,
#' it is possible that these SNPs could be recovered after an examination of the results. Be
#' aware that testing of allele frequencies might fail for small datasets.
#'
#' @export
#' @param SNPobj List with SNPs data from plink or \code{VCF-class}.
#' @param checkAlleleFreqs Should allele frequencies be check (Default: TRUE)
#' @return List containing the SNPs prepared for \code{scoreInvHap}
#' \itemize{
#' \item{genos: Object with genotype data ready for scoreInvHap}
#' \item{wrongAlleles: Character vector with the SNPs discarded due to having alleles different to reference}
#' \item{wrongFreqs: Character vector with the SNPs discarded due to having allele frequencies different to reference}
#' }
#' @examples
#'
#' ## Run method
#' if(require(VariantAnnotation)){
#'     vcf <- readVcf(system.file("extdata", "example.vcf", package = "scoreInvHap"), "hg19")
#'     resList <- checkSNPs(vcf)
#'     resList
#' }
checkSNPs <- function(SNPobj, checkAlleleFreqs = TRUE){

    info <- scoreInvHap::info

    ## Check VCF
    if (is(SNPobj, "VCF")){
        dupSNPs <- rownames(SNPobj)[duplicated(rownames(SNPobj))]
        SNPobj <- SNPobj[!rownames(SNPobj) %in% dupSNPs, ]

        ## Remove variants with length != 1
        SNPobj <- SNPobj[width(SNPobj) == 1, ]

        ## Filter SNPs with bad imputation quality
        if ("R2" %in% colnames(VariantAnnotation::info(SNPobj))){
            SNPobj <- SNPobj[VariantAnnotation::info(SNPobj)$R2 > 0.4, ]
        }

        ranges <- SummarizedExperiment::rowRanges(SNPobj)
        rownames(SNPobj) <- paste(seqnames(ranges), start(ranges), sep = ":")

        ## Select SNPs in scoreInvHap
        SNPobj <- SNPobj[rownames(SNPobj) %in% rownames(info), ]

        ## Check alleles
        ranges <- SummarizedExperiment::rowRanges(SNPobj)
        ref <- as.character(ranges$REF)
        alt <- as.character(unlist(ranges$ALT))

        info <- info[rownames(SNPobj),]
        badMask <- (ref != info$Ref & ref != info$Alt) | (alt != info$Ref & alt != info$Alt)
        badSNPs <- NULL

        ## Remove SNPs with wrong alleles
        if (sum(badMask) > 0){
            badSNPs <- rownames(SNPobj)[badMask]
            SNPobj <- SNPobj[!badMask, ]
            info <- info[!badMask, ]
        }

        ## Check allele Freqs
        if(checkAlleleFreqs){

            ranges <- SummarizedExperiment::rowRanges(SNPobj)
            ref <- as.character(ranges$REF)

            stats <- VariantAnnotation::snpSummary(SNPobj)
            freq <- ifelse(ref == info$Ref, stats$a1Freq, stats$a2Freq)
            freqMask <- abs(freq - info$Freq) > 0.2

            badFreq <- NULL
            if (sum(freqMask) > 0){
                badFreq <- rownames(SNPobj)[freqMask]
                SNPobj <- SNPobj[!freqMask, ]
                info <- info[!freqMask, ]
            }
        }

        ## Change ids
        rownames(SNPobj) <- info[rownames(SNPobj), "id"]
        return(list(genos = SNPobj, wrongAlleles = badSNPs,  wrongFreqs = badFreq))
    }

    ## Check plink
    if (!all(c("map", "genotypes") %in% names(SNPobj))){
        stop("SNPobj must contain map and genotypes elements")
    }

    map <- SNPobj$map
    geno <- SNPobj$geno

    if (!all(c("allele.1", "allele.2") %in% colnames(map))){
        stop("map of SNPobj must contain columns allele.1 and allele.2")
    }

    map <- map[!is.na(map$allele.1), ]

    ## Change name to position
    if (sum(duplicated(paste(map$chromosome, map$position, sep = ":")))){
        stop("There are different SNPs located in the same position. Please, include only one variant per site before proceeding.")
    }
    colnames(geno) <- rownames(map) <- paste(map$chromosome, map$position, sep = ":")

    ## Select SNPs in scoreInvHap
    selSNPs <- rownames(map)[rownames(map) %in% rownames(info)]
    map <- map[selSNPs, ]
    geno <- geno[, selSNPs]

    ## Check alleles
    info <- info[rownames(map),]
    badMask <- (map$allele.1 != info$Ref & map$allele.1 != info$Alt) | (map$allele.2 != info$Ref & map$allele.2 != info$Alt)
    badSNPs <- NULL

    ## Remove SNPs with wrong alleles
    if (sum(badMask) > 0){
        badSNPs <- rownames(map)[badMask]
        map <- map[!badMask, ]
        geno <- geno[, !badMask]
        info <- info[!badMask, ]
    }


    ## Check allele Freqs
    if(checkAlleleFreqs){

        stats <- snpStats::col.summary(geno)
        freq <- ifelse(map$allele.1 == info$Ref, stats$RAF, 1 - stats$RAF)
        freqMask <- abs(freq - info$Freq) > 0.2

        badFreq <- NULL
        if (sum(freqMask) > 0){
            badFreq <- rownames(map)[freqMask]
            map <- map[!freqMask, ]
            geno <- geno[, !freqMask]
            info <- info[!freqMask, ]

        }
    }

    ## Change ids
    colnames(geno) <- rownames(map) <- info[rownames(map), "id"]

    SNPobj$map <- map
    SNPobj$genotypes <- geno

    return(list(genos = SNPobj, wrongAlleles = badSNPs,  wrongFreqs = badFreq))
}
