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

test_writeVcf_tags <- function()
{
    fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
    dest <- tempfile()
    vcf1 <- readVcf(fl, "hg19")
    hd1 <- exptData(vcf1)$header
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    hd2 <- exptData(vcf2)$header
    checkTrue(rownames(meta(hd1)) %in% rownames(meta(hd2))) 
    checkIdentical(names(geno(vcf1)), names(geno(vcf2))) 
    checkIdentical(colnames(mcols(info(vcf1))), colnames(mcols(info(vcf2))))
}
 
test_writeVcf_flatgeno <- function()
{
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    dest <- tempfile()
    vcf1 <- readVcf(fl, "hg19")
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
}
