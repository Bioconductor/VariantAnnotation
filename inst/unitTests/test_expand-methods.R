expand <- VariantAnnotation::expand 

test_expand_info_geno <- function()
{
    fl <- system.file("unitTests", "cases", "expand.vcf", 
        package="VariantAnnotation")
    vcf <- suppressWarnings(readVcf(fl, "hg19"))
    exp <- expand(vcf)
    names(mcols(rowData(exp)))
    cnms <- c("paramRangeID", "REF", "ALT", "QUAL", "FILTER")
    checkIdentical(cnms, names(mcols(rowData(exp))))
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

test_expand_structural <- function()
{
    fl <- system.file("unitTests", "cases", "mixedStructural.vcf", 
                      package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    exp <- expand(vcf)
    checkTrue(nrow(exp) == 8L)
}

test_expand_multiple_info <- function()
{
    fl <- system.file("unitTests", "cases", "multiple_INFO_fields.vcf", 
                      package="VariantAnnotation")
    vcf <- readVcf(fl, "GRCh37")
    exp <- expand(vcf)
    checkTrue(nrow(exp) == 3L)

    ## single row, subset of INFO fields selected 
    exp <- expand(vcf[1,])
    checkTrue(nrow(exp) == 2L)

    ## single row, single (only) INFO field selected
    vcf <- readVcf(fl, "GRCh37", param=ScanVcfParam(info="AF"))
    exp <- expand(vcf[1,])
    checkTrue(nrow(exp) == 2L)
}
