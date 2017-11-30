#'#################################################################################
#'#################################################################################
#' Create scoreInvHap objects
#'#################################################################################
#'#################################################################################

#'#################################################################################
# Load data ####
#'#################################################################################

#' This section contains the code used to load SNP data from VCF using functions
#' defined in extraFunctions.R.

# Load packages and scripts  ####
#'#################################################################################
library(GenomicRanges)
library(VariantAnnotation)
library(parallel)

## Load functions from extraFunctions file
source("./inst/scripts/extraFunctions.R")

# Load ancestry
load("/SYNCRW10125/DATASETS/STUDY/1000GENOME/Samples_Pop1GK.Rdata")

# Load genotypes
load("/DATA/DATAW5/Carlos/Inversions/SNPfierTesting/invFESTinvClustcomplex.Rdata")
load("/home/cruiz/InversionSequencing/SandersROIs/invClust/invClustEURNoSeqDups.Rdata")

# Create ranges
load("/home/cruiz/InversionSequencing/SandersROIs/ROIsGR.Rdata")

ranges <- ROIsGR[c("ROIno.8.3", "ROIno.17.16")]

ranges <- c(granges(ranges),
            GRanges(c("7:54290974-54386821", "X:72215927-72306774")))
names(ranges)[3:4] <- names(invFestGenotypescomplex)

# Load SNPs  ####
#'#################################################################################

selinvs <- names(ranges)

## Select samples used to create references
EUR <- rownames(samp_pop)[samp_pop$superpop == "EUR"]


# Get region SNPs

# Set vcffile parameter to path to vcf file with the genotypes
inversionSNPs <- list(getVCFmatrix(ranges[1], EUR, minmaf = 0.05),
                      getVCFmatrix(ranges[2], EUR, minmaf = 0.05),
                      getVCFmatrix(ranges[3], EUR, minmaf = 0.01),
                      getVCFmatrix(ranges[4], EUR, minmaf = 0.01,
                                   vcffile = "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"))
names(inversionSNPs) <- selinvs


#'#################################################################################
# Generate objects ####
#'#################################################################################

#' This section contains the code used to generate the references from genotypes
#' data and inversion classification.

# Objects decription
## - invsGeno/invFestGenotypescomplex: list with the inversion status of the samples
##                                     (each element is a factor with inversion calling)
## - inversionSNPs: list with SNP data. Each element is a list with two elements:
##      - genotypes: SnpMatrix with samples genotypes
##      - map: data.frame with SNPs annotation
## - selinvs: names of the inversions (used to run lapply)

# Compute SNPs R2 with the inversion
SNPsR2 <- c(
    lapply(selinvs[1:2], function(x){
        inv <- as.numeric(invsGeno[[x]])
        inv <- new("SnpMatrix", as.matrix(inv))
        ld(inv, inversionSNPs[[x]]$genotypes, stats=c("R.squared"))[1,]
    }), list({
        haplos <- c("Ia", "Ib", "Na", "Nb")
        R2s <- vapply(haplos, function(hap) {
            pats <- gregexpr(hap, invFestGenotypescomplex[["HsInv0286"]])
            lenvec <- lengths(pats)
            lenvec[lenvec == 1 & sapply(pats, function(y) all(y == -1))] <- 0
            inv <- new("SnpMatrix", as.matrix(lenvec + 1))
            ld(inv, inversionSNPs[["HsInv0286"]]$genotypes, stats=c("R.squared"))[1,]
        }, numeric(ncol(inversionSNPs[["HsInv0286"]]$genotypes)))
        res <- matrixStats::rowMaxs(R2s)
        names(res) <- rownames(R2s)
        res
    }, {
        haplos <- c("I", "Nc", "Na", "Nb")
        R2s <- vapply(haplos, function(hap) {
            pats <- gregexpr(hap, invFestGenotypescomplex[["HsInv0396"]])
            lenvec <- lengths(pats)
            lenvec[lenvec == 1 & sapply(pats, function(y) all(y == -1))] <- 0
            inv <- new("SnpMatrix", as.matrix(lenvec + 1))
            ld(inv, inversionSNPs[["HsInv0396"]]$genotypes, stats=c("R.squared"))[1,]
        }, numeric(ncol(inversionSNPs[["HsInv0396"]]$genotypes)))
        res <- matrixStats::rowMaxs(R2s)
        names(res) <- rownames(R2s)
        res
    }))
names(SNPsR2) <- selinvs

# Compute allele frequency for the different haplotypes
Refs <- c(lapply(selinvs[1:2], function(x)
  computeRefs(invsGeno[[x]], inversionSNPs[[x]]$genotypes,
              inversionSNPs[[x]]$map)),
  lapply(selinvs[3:4], function(x)
    computeRefs(invFestGenotypescomplex[[x]], inversionSNPs[[x]]$genotypes,
                inversionSNPs[[x]]$map)))
names(Refs) <- selinvs

## Create a character vector with the heterozygote form in 1000Genomes to ensure
## compatibility with other arrays
hetRefs <- lapply(selinvs, function(x){
  map <- inversionSNPs[[x]]$map
  rownames(map) <- map[,1]
  genos <- inversionSNPs[[x]]$genotypes
  alleletable <- getAlleleTable(map)
  res <- alleletable[, 2]
  names(res) <- rownames(alleletable)
  res
})
names(hetRefs) <- selinvs

save(hetRefs, Refs, SNPsR2,
     file = "./newSNPfierdata.Rdata")
