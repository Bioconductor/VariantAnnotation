test_snpSummary_empty <- function()
{
    cnames <- c("g00", "g01", "g11", "a0Freq", "a1Freq", 
                "HWEzscore", "HWEpvalue")
    ## no genotype
    fl <- system.file("unitTests", "cases", 
              "FORMAT_header_no_SAMPLEs.vcf", 
              package="VariantAnnotation")
    res <- suppressWarnings(snpSummary(readVcf(fl, "hg19")))
    checkTrue(all(dim(res) == c(0L, 7L)))
    checkEquals(names(res), cnames)

    ## structural
    fl <- system.file("extdata", "structural.vcf", 
              package="VariantAnnotation")
    res <- suppressWarnings(snpSummary(readVcf(fl, "hg19")))
    checkTrue(all(dim(res) == c(0L, 7L)))
    checkEquals(names(res), cnames)
}

test_snpSummary_output <- function()
{
    fl <- system.file("extdata", "ex2.vcf", 
              package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    res <- snpSummary(vcf)

    checkTrue(is.integer(res$g00))
    checkTrue(is.integer(res$g01))
    checkTrue(is.integer(res$g11))

    checkEquals(nrow(res), nrow(vcf))
    checkEquals(rownames(res), rownames(vcf))

    checkEquals(as.integer(res[1,c("g00", "g01", "g11")]), c(1L, 1L, 1L))
    checkEquals(as.integer(res[2,c("g00", "g01", "g11")]), c(2L, 1L, 0L))
    checkEquals(round(res[["a0Freq"]], 2), c(0.50, 0.83, NA, NA, NA))
    checkEquals(round(res[["a1Freq"]], 2), c(0.50, 0.17, NA, NA, NA))
}

