
test_expand_info_geno <- function()
{
    fl <- system.file("unitTests", "cases", "expand.vcf", 
        package="VariantAnnotation")
    vcf <- suppressWarnings(readVcf(fl, "hg19"))
    exp <- expand(vcf)
    checkTrue(nrow(exp) == 22L)
    checkIdentical(info(exp)$DP[5:6], c(10L, 10L))
    checkIdentical(info(exp)$AF[5:6], c(0.333, 0.667))
    checkIdentical(geno(exp)$AAP[6:11], c(2L, rep(1L, 4L), 3L))
}

test_expand_structural <- function()
{
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    exp <- expand(vcf)
    checkTrue(nrow(vcf) == nrow(exp))
    checkTrue(ncol(vcf) == ncol(exp))
    checkTrue(is.character(alt(exp)))
}
