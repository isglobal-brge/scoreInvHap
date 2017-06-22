#'#################################################################################
#'#################################################################################
#' Create scoreInvHap objects
#'#################################################################################
#'#################################################################################

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

# Compute SNPs R2 with the inversion
library(rms)
SNPsR2 <- c(
  lapply(selinvs[1:2], function(x){
  genos <- as(inversionSNPs[[x]]$genotypes, "character")
  ld <- apply(genos, 2, function(snps) {
    df <- data.frame(inv = invsGeno[[x]], snp = snps)
    lrm(inv ~ snp, df, maxit = 100, tol = 1e-29)$stats["R2"]
  })
  ld
}),
lapply(selinvs[3:4], function(x){
  genos <- as(inversionSNPs[[x]]$genotypes, "character")
  ld <- apply(genos, 2, function(snps) {
    df <- data.frame(inv = invFestGenotypescomplex[[x]], snp = snps)
    lrm(inv ~ snp, df, maxit = 100, tol = 1e-29)$stats["R2"]
  })
  ld
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
