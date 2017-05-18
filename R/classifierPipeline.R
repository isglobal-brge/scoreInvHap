#' Classify inversions using identity
#'
#' This is the main function of `snpfier` package. This function accepts
#' SNPs data in a plink or a VCF format and compute the inversion prediction.
#'
#' @export
#' @param SNPlist List with SNPs data. It should contain genotypes (a \code{SNPmatrix}) and map (a data.frame
#' with the annotation)
#' @param SNPsR2 Vector with the R2 of the SNPs of the region
#' @param hetRefs Vector with the heterozygote form of the SNP in the inversion
#' @param Refs List with the allele frequencies in the references
#' @param R2 Vector with the R2 between the SNPs and the inversion status
#' @param imputed Logical. If TRUE, scores are computed using posterior probabilities.
#' If FALSE, scores are computed using best guess. Only applied when SNPlist is a VCF.
#' @param mc.cores Numeric with the number of cores used in the computation
#' @param verbose Should message be shown?
#' @return A \code{SNPfieRes} object
#' @examples
#' if(require(VariantAnnotation)){
#'   vcf <- readVcf(system.file("extdata", "example.vcf", package = "snpfier"), "hg19")
#'   res <- classifierPipeline(vcf, SNPsR2$HsInv0286, hetRefs = hetRefs$HsInv0286,
#'   Refs$HsInv0286, mc.cores = 1)
#' }
classifierPipeline <- function(SNPlist, SNPsR2, hetRefs, Refs, R2 = 0,
                               imputed = FALSE,
                               mc.cores = 1, verbose = FALSE){
  if (is(SNPlist, "VCF")){
    dupSNPs <- rownames(SNPlist)[duplicated(rownames(SNPlist))]
    SNPlist <- SNPlist[!rownames(SNPlist) %in% dupSNPs, ]

    if (imputed){

      ## Select SNPs with a R2 equal or higher than the threshold
      SNPsR2 <- SNPsR2[SNPsR2 >= R2]

      ## Filter objects to only those included in the references
      commonSNPs <- Reduce(intersect, list(names(SNPsR2), names(hetRefs), names(Refs), rownames(SNPlist)))

      if (!length(commonSNPs)){
        stop("There are no common SNPs between the SNP object and the reference.")
      }

      SNPlist <- SNPlist[commonSNPs, ]

      ## Filter SNPs with bad imputation quality
      SNPlist <- SNPlist[info(SNPlist)$R2 > 0.4, ]

      # Create alleleTable
      map <- prepareMap(SNPlist)
      alleletable <- getAlleleTable(map)
      alleletable <- correctAlleleTable(alleletable = alleletable, hetRefs = hetRefs, map = map)

      genos <- geno(SNPlist)$GP

      haploid <- ifelse(dim(genos)[3] == 2, TRUE, FALSE)

      # Order alleles in Refs as in our VCF
      refs <- adaptRefs(Refs = Refs, alleletable = alleletable, haploid = haploid)


      ## Compute score
      classifScore <- classifSNPsImpute(genos, R2 = SNPsR2, refs = refs)
      inv <- getInvStatus(scores = classifScore$scores)
      res <- new("SNPfieRes", classification = inv$class,
                 scores = classifScore$scores, numSNPs = classifScore$numSNPs,
                 certainty = inv$certainty)
      return(res)
    } else {
      vcf <- SNPlist
      SNPlist <- VariantAnnotation::genotypeToSnpMatrix(vcf)

      SNPlist$map$position <- start(rowRanges(vcf))
      SNPlist$map$chromosome <- as.character(seqnames(rowRanges(vcf)))
      SNPlist$genotypes <- SNPlist$genotypes[, !SNPlist$map$ignore]
      SNPlist$map <- SNPlist$map[!SNPlist$map$ignore, ]
      SNPlist$map$allele.2 <- unlist(SNPlist$map$allele.2)
      rownames(SNPlist$map) <- SNPlist$map$snp.names
    }
  }

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
  commonSNPs <- Reduce(intersect, list(names(SNPsR2), names(hetRefs), names(Refs), rownames(map)))

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
  geno <- SNPlist$genotypes[, rownames(alleletable)]
  genos <- getGenotypesTable(geno = geno, allele = alleletable)

  if (verbose){
    message("Computing scores")
  }
  classifScore <- classifSNPs(genos = genos, R2 = SNPsR2, refs = Refs,
                              mc.cores = mc.cores)
  inv <- getInvStatus(scores = classifScore$scores)
  res <- new("SNPfieRes", classification = inv$class,
             scores = classifScore$scores, numSNPs = classifScore$numSNPs,
             certainty = inv$certainty)
  return(res)
}


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

prepareMap <- function(vcf){
  map <- GenomicRanges::mcols(rowRanges(vcf))
  rownames(map) <- rownames(vcf)
  cnmap <- colnames(map)
  cnmap[cnmap == "REF"] <- "allele.1"
  cnmap[cnmap == "ALT"] <- "allele.2"
  colnames(map) <- cnmap
  map$allele.2 <- unlist(map$allele.2)
  map
}
