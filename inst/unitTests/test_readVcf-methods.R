f1 <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
f2 <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")

test_readVcf <- function()
{
    vcf <- readVcf(f1)
    checkIdentical(names(assays(vcf)), c("GT", "GQ", "DP"))
    checkEquals(width(rowData(vcf)), width(values(rowData(vcf))[["REF"]])) 
    checkEquals(length(unlist(values(rowData(vcf))[["ALT"]])), 18) 
    checkEquals(scanVcf(f1)[[1]]$POS, start(rowData(vcf))) 
}
