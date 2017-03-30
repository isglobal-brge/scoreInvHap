#' Classify inversions using identity
#'
#' @export
#' @param SNPlist List with SNPs data. It should contain genotypes (a SNPmatrix) and map (a data.frame
#' with the annotation)
#' @param SNPsR2 Vector with the R2 of the SNPs of the region
#' @param hetRefs Vector with the heterozygote form of the SNP in the inversion
#' @param Refs List with the allele frequencies in the references
#' @param R2 Vector with the R2 between the SNPs and the inversion status
#' @param mc.cores Numeric with the number of cores used in the computation
#' @param verbose Should message be shown?
#' @return A list with the following elements:
#' \itemize{
#' \item{class: Factor with the classification of the individuals}
#' \item{scores: Matrix with the simmilarity scores of the individuals}
#' }
classifierPipeline <- function(SNPlist, SNPsR2, hetRefs, Refs, R2 = 0.3,
                               alfreq, genofreq,
                               mc.cores = 1, verbose = FALSE){
  if (!all(c("map", "genotypes") %in% names(SNPlist))){
    stop("SNPlist must contain map and genotypes elements")
  }

  ## Prevent errors with forking in windows
  if( .Platform$OS.type == "windows" ){
    mc.cores <- 1
  }

  map <- SNPlist$map

  if (!all(c("allele.1", "allele.2") %in% colnames(map))){
    stop("map of SNPlist must contain columns allele.1 and allele.2")
  }

  map <- map[!is.na(map$allele.1), ]


  ## Select SNPs with a R2 equal or higher than the threshold
  SNPsR2 <- SNPsR2[SNPsR2 >= R2]

  ## Filter objects to only those included in the references
  commonSNPs <- Reduce(intersect, list(names(SNPsR2), names(hetRefs), names(refs), rownames(map), names(alfreq)))

  if (!length(commonSNPs)){
    stop("There are no common SNPs between the SNP object and the reference.")
  }

  map <- map[commonSNPs, ]

  if (verbose){
    message("Computing allele table")
  }
  alleletable <- getAlleleTable(map = map)
  alleletable <- correctAlleleTable(alleletable = alleletable, hetRefs = hetRefs, map = map)

  if (verbose){
    message("Computing genotype table")
  }
  geno <- SNPlist$genotypes[, rownames(map)]
  genos <- getGenotypesTable(geno = geno, allele = alleletable)

  if (verbose){
    message("Computing scores")
  }
  classifScore <- classifSNPspar(genos = genos, R2 = SNPsR2, refs = Refs,
                                alfreq = alfreq, genofreq = genofreq,
                                mc.cores = mc.cores)
  inv <- getInvStatus(scores = classifScore$scores)
  res <- new("SNPfieRes", classification = inv$class, certainty = classifScore$probs,
             scores = classifScore$scores, numSNPs = classifScore$numSNPs)
  res
}
