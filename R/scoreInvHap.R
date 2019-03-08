#' Classify inversions using identity
#'
#' This is the main function of `scoreInvHap` package. This function accepts
#' SNPs data in a plink or a VCF format and compute the inversion prediction. The list
#' of available inversions is included in a GenomicRanges called `inversionGR`.
#'
#' @export
#' @param SNPlist List with SNPs data from plink or \code{VCF-class}.
#' @param inv Character with the name of the inversion to genotype. The available
#' inversions are included in a table in the main vignette.
#' @param SNPsR2 Vector with the R2 of the SNPs of the region
#' @param hetRefs Vector with the heterozygote form of the SNP in the inversion
#' @param Refs List with the allele frequencies in the references
#' @param R2 Vector with the R2 between the SNPs and the inversion status
#' @param probs Logical. If TRUE, scores are computed using posterior probabilities.
#' If FALSE, scores are computed using best guess. Only applied when SNPlist is a VCF.
#' @param BPPARAM A \code{BiocParallelParam} instance. Used to parallelize computation
#' @param verbose Should message be shown?
#' @return A \code{scoreInvHap} object
#' @examples
#'
#' # See list of inversions
#' data(inversionGR)
#' inversionGR
#'
#' ## Run method
#' if(require(VariantAnnotation)){
#'     vcf <- readVcf(system.file("extdata", "example.vcf", package = "scoreInvHap"), "hg19")
#'     res <- scoreInvHap(vcf, inv = "inv7_005")
#' }
#'
#'
scoreInvHap <- function(SNPlist, inv = NULL, SNPsR2, hetRefs, Refs, R2 = 0,
                        probs = FALSE,
                        BPPARAM = BiocParallel::SerialParam(), verbose = FALSE){

    if (is.character(inv)){
        SNPsR2 <- scoreInvHap::SNPsR2[[inv]]
        hetRefs <- scoreInvHap::hetRefs[[inv]]
        Refs <- scoreInvHap::Refs[[inv]]
    }

    if (is(SNPlist, "VCF")){
        if (imputed){

            ## Select SNPs with a R2 equal or higher than the threshold
            SNPsR2 <- SNPsR2[SNPsR2 >= R2]

            ## Filter objects to only those included in the references
            commonSNPs <- Reduce(intersect, list(names(SNPsR2), names(hetRefs),
                                                 names(Refs), rownames(SNPlist)))

            if (!length(commonSNPs)){
                stop("There are no common SNPs between the SNP object and the reference.")
            }

            SNPlist <- SNPlist[commonSNPs, ]

            # Create alleleTable
            map <- prepareMap(SNPlist)
            alleletable <- getAlleleTable(map)

            genos <- geno(SNPlist)$GP

            haploid <- ifelse(dim(genos)[3] == 2, TRUE, FALSE)

            # Order alleles in Refs as in our VCF
            refs <- adaptRefs(Refs = Refs, alleletable = alleletable,
                              haploid = haploid)


            ## Compute score
            classifScore <- classifSNPsImpute(genos, R2 = SNPsR2, refs = refs,
                                              BPPARAM = BPPARAM)
            inv <- getInvStatus(scores = classifScore$scores)
            res <- new("scoreInvHapRes", classification = inv$class,
                       scores = classifScore$scores, numSNPs = classifScore$numSNPs,
                       certainty = inv$certainty)
            return(res)
        } else {
            vcf <- SNPlist
            SNPlist <- VariantAnnotation::genotypeToSnpMatrix(vcf)

            SNPlist$map$position <- GenomicRanges::start(SummarizedExperiment::rowRanges(vcf))
            SNPlist$map$chromosome <- as.character(GenomicRanges::seqnames(SummarizedExperiment::rowRanges(vcf)))
            SNPlist$genotypes <- SNPlist$genotypes[, !SNPlist$map$ignore]
            SNPlist$map <- SNPlist$map[!SNPlist$map$ignore, ]
            SNPlist$map$allele.2 <- unlist(SNPlist$map$allele.2)
            rownames(SNPlist$map) <- SNPlist$map$snp.names
        }
    }

    if (!all(c("map", "genotypes") %in% names(SNPlist))){
        stop("SNPlist must contain map and genotypes elements")
    }

    map <- SNPlist$map

    if (!all(c("allele.1", "allele.2") %in% colnames(map))){
        stop("map of SNPlist must contain columns allele.1 and allele.2")
    }

    map <- map[!is.na(map$allele.1), ]


    ## Select SNPs with a R2 equal or higher than the threshold
    SNPsR2 <- SNPsR2[SNPsR2 >= R2]

    ## Filter objects to only those included in the references
    commonSNPs <- Reduce(intersect, list(names(SNPsR2), names(hetRefs), names(Refs), rownames(map)))

    if (!length(commonSNPs)){
        stop("There are no common SNPs between the SNP object and the reference.")
    }

    map <- map[commonSNPs, ]

    if (verbose){
        message("Computing allele table")
    }
    alleletable <- getAlleleTable(map = map)

    if (verbose){
        message("Computing genotype table")
    }
    geno <- SNPlist$genotypes[, rownames(alleletable)]

    if (verbose){
        message("Computing scores")
    }
    classifScore <- classifSNPs(genos = geno, R2 = SNPsR2, refs = Refs,
                                alleletable = alleletable,
                                BPPARAM = BPPARAM)
    inv <- getInvStatus(scores = classifScore$scores)
    res <- new("scoreInvHapRes", classification = inv$class,
               scores = classifScore$scores, numSNPs = classifScore$numSNPs,
               certainty = inv$certainty)
    return(res)
}

#' Adapt references to imputed data
#'
#' @description Internal
#'
#' @param Refs List with the allele frequencies
#' @param alleletable Data.frame with the alleles per SNP (from getAlleleTable)
#' @param haploid Logical. If TRUE, modify references for haploid samples
#' @return List with the same values than Refs but adapted to imputation data
adaptRefs <- function(Refs, alleletable, haploid = FALSE){
    rb <- lapply(rownames(alleletable), function(snp) {
        als <- unlist(alleletable[snp, 1:3])
        rb <- Refs[[snp]]
        if (ncol(rb) == 1){
            rb <- cbind(rb, 0, 0)
            colnames(rb)[2:3] <- als[!als %in% colnames(rb)]
        }
        if (ncol(rb) == 2){
            rb <- cbind(rb, 0)
            colnames(rb)[3] <- als[!als %in% colnames(rb)]
        }
        rb <- rb[, als]
        if (haploid){
            rb <- rb[, -2]
        }
        rb
    })
    names(rb) <- rownames(alleletable)
    rb
}

#' Modify feature data from VCF
#'
#' @description Internal. Modify feature data from VCF to comply with
#' scoreInvHap requirements.
#'
#' @param vcf \code{VCF} object
#' @return Data.frame with the feature data
prepareMap <- function(vcf){
    map <- GenomicRanges::mcols(SummarizedExperiment::rowRanges(vcf))
    rownames(map) <- rownames(vcf)
    cnmap <- colnames(map)
    cnmap[cnmap == "REF"] <- "allele.1"
    cnmap[cnmap == "ALT"] <- "allele.2"
    colnames(map) <- cnmap
    map$allele.2 <- unlist(map$allele.2)
    map
}
