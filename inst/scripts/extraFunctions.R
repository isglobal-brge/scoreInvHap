#'#################################################################################
#'#################################################################################
#' Functions required to create scoreInvHap references
#'#################################################################################
#'#################################################################################

#' Get SNP matrix from VCF file
#' @param range GenomicRanges with the range of the whole regions
#' @param samples Character with the samples to be selected in the VCF
#' @param snps.names Character with the snps to be selected in the VCF
#' @param minmaf Numeric with the minimum maf required for a position to be in the final snpMatrix
getVCFmatrix <- function(range, samples = NULL, snps.names = NULL, minmaf = 0.1, Remove.Granges = NULL, ...){
    vcf <- loadVCFrange(range, samples, ...)
    vcf <- filterVCF(vcf, snps.names = snps.names, Remove.Granges = Remove.Granges)

    # Convert haploid genotypes to diploid (for males in chrX)
    geno(vcf)$GT <- sapply(1:ncol(geno(vcf)$GT), function(x) {
        col <-  geno(vcf)$GT[, x]
        if(nchar(col[1]) == 1){
            col <- paste0(col, "|", col)
        }
        col
    })
    snpsVCF <- genotypeToSnpMatrix(vcf)
    ## Conversion from VCF to SNP matrix produces some SNPs to be NA (multiallelic or bigger than 1)
    snpsVCF$map$position <- start(rowRanges(vcf))
    snpsVCF$map$chromosome <- rep(as.character(seqnames(range)), nrow(snpsVCF$map))
    snpsVCF$genotypes <- snpsVCF$genotypes[, !snpsVCF$map$ignore]
    snpsVCF$map <- snpsVCF$map[!snpsVCF$map$ignore, ]
    snpsVCF$map$allele.2 <- unlist(snpsVCF$map$allele.2)
    sums <- col.summary(snpsVCF$genotypes)
    snpsVCF$map <- snpsVCF$map[!is.na(sums$MAF), ]
    snpsVCF$genotypes <- snpsVCF$genotypes[, !is.na(sums$MAF)]

    sums <- col.summary(snpsVCF$genotypes)
    inds <- sums$MAF > minmaf
    snpsVCF$genotypes <- snpsVCF$genotypes[, inds]
    snpsVCF$map <- snpsVCF$map[inds, ]
    rownames(snpsVCF$map) <- snpsVCF$map$snp.names
    snpsVCF
}


#' @param range GenomicRanges with the range of the VCF to be loaded
#' @param samples Character vector with the samples to be loaded
#' @return A collapsedVCF object
loadVCFrange <- function(range, samples = NULL,
                         vcffile = paste0("ALL.chr", as.character(seqnames(range)), ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")){
    vcfsamples <- samples(scanVcfHeader(vcffile))
    if (!is.null(samples)){
        samples <- vcfsamples[vcfsamples %in% samples]
    } else{
        samples <- vcfsamples
    }
    param <- ScanVcfParam(samples = samples, which = range)
    vcf <- readVcf(vcffile, "hg19", param)
    vcf
}

#' @param vcf A collapsedVCF object
#' @param snps.names Character vector with the SNPs that are selected
#' @return A collapsedVCF object without multiallelic variants
filterVCF <- function(vcf, snps.names = NULL, Remove.Granges = NULL){
    if (!is.null(snps.names)){
        vcf <- vcf[rownames(vcf) %in% snps.names, ]
    }
    if (!is.null(Remove.Granges)){
        hits <- findOverlaps(rowRanges(vcf), Remove.Granges)
        if (length(hits)){
            vcf <- vcf[-from(hits),  ]
        }
    }

    ### Remove multiallelic variants
    vcf_filt <- vcf[lengths(alt(vcf)) == 1]
    vcf_filt
}


computeRefs <- function(haplos, genos, map){
    rownames(map) <- map$snp.names
    map <- map[colnames(genos), ]
    alleletable <- getAlleleTable(map)
    geno <- getGenotypesTable(genos, alleletable)

    Refs <- lapply(colnames(geno), function(x) {
        t <- as(table(haplos, geno[, x]), "matrix")
        if ("NN" %in% colnames(t)){
            t <- t[-which(colnames(t) == "NN")]
        }
        t <- t/rowSums(t)
        t
    })
    names(Refs) <- colnames(geno)
    Refs
}



#' @param map data.frame with the annotation of the SNPs (from plink format)
getAlleleTable <- function(map){
    AB <- matrix(c(as.character(map$allele.1), as.character(map$allele.2)), ncol = 2)
    allele <- data.frame(AA = paste0(map$allele.1, map$allele.1),
                         AB = apply(AB, 1, function(x) paste(sort(x), collapse = "")),
                         BB = paste0(map$allele.2, map$allele.2), stringsAsFactors = FALSE)
    rownames(allele) <- rownames(map)
    allele$Fail <- rep("NN", nrow(allele))
    allele
}

#' @param geno SnpMatrix (from plink format)
#' @param allele data.frame with the alleles per SNP (from getAlleleTable)
getGenotypesTable <- function(geno, allele){
    names <- rownames(geno)

    ## Convert matrix of genotypes to numeric (0, 1, 2). Add 1 to match column indeces of allele table.
    geno <- as(geno, "numeric") + 1
    geno[is.na(geno)] <- 4
    geno <- vapply(1:ncol(geno),function(x) unlist(allele[x, geno[,x], drop = TRUE]), character(nrow(geno)))
    colnames(geno) <- rownames(allele)
    rownames(geno) <- names
    geno
}
