test_writeVcf_connection_increment <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf1 <- readVcf(fl, "hg19")

    outfl <- tempfile()
    con <- file(outfl, open="a")
    writeVcf(vcf1[1:2,], con)
    writeVcf(vcf1[-(1:2),], con)
    close(con)
    vcf2 <- readVcf(outfl, "hg19")

    checkIdentical(dim(vcf1), dim(vcf2))
}
