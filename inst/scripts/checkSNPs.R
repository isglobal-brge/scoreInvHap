#'#################################################################################
#' Check SNPs before running scoreInvHap
#'#################################################################################

## Create List with SNPs ids in scoreInvHap (R)
library(scoreInvHap)
ids <- names(Reduce(c, SNPsR2))
write.table(ids, file = "scoreInvHap_SNPs.txt", quote = FALSE, row.names = FALSE,
            col.names = FALSE)

## Create variant list file from 1000 Genomes
vcftools --gzvcf ~/PublicData/STUDY/1000GENOME/VCF/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --snps scoreInvHap_SNPs.txt --recode --recode-INFO EUR_AF --out scoreInvHap.sites

## Create list positions of selected SNPs
grep -vE "^#" scoreInvHap.sites.recode.vcf | awk '{print $1, $2, $2, $3}' > selSNPsPos.tab

# Filter VCF
vcftools --gzvcf test.vcf --positions selSNPsPos.tab --recode --out plink

# Filter plink
plink -bfile newAnnot --extract range selSNPsPos.tab --make-bed --out plink


## Create data.frame with scoreInvHap SNPs information (R)
info <- read.table("scoreInvHap.sites.recode.vcf", as.is = TRUE)
info <- info[, -c(6,7)]
colnames(info) <- c("chromosome", "position", "id", "Ref", "Alt", "Freq")
info$Freq <- as.numeric(gsub("EUR_AF=", "", info$Freq))
rownames(info) <- paste(info$chromosome, info$position, sep = ":")


