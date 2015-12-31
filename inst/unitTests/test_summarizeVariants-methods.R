grl <- GRangesList(
    A=GRanges("20", IRanges(14370, 17331)),
    B=GRanges("20", IRanges(17330, width=1)),
    C=GRanges("20", IRanges(1110696, width=1)), 
    D=GRanges("20", IRanges(1110696, width=2)))

test_summarizeVariants <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    sv <- summarizeVariants(grl, vcf, findOverlaps)
    checkIdentical(ncol(vcf), ncol(sv))
    checkIdentical(length(grl), nrow(sv))

    target <- matrix(c(0, 0, 1, 1, 2, 1, 1, 1, 1, 0, 1, 1), ncol=3)
    current <- unname(assays(sv)$counts)
    checkIdentical(current, target) 
}

