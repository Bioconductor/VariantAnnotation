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

test_writevcf_geno <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    dest <- tempfile()

    ## empty
    vcf1 <- readVcf(fl, "hg19", param=ScanVcfParam(geno=NA))
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkTrue(length(geno(vcf2)$GT) == 0L)

    ## matrix 
    vcf1 <- readVcf(fl, "hg19", param=ScanVcfParam(geno="GT"))
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT)

    param=ScanVcfParam(geno="GT", samples="NA00003")
    vcf1 <- readVcf(fl, "hg19", param=param)
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT)

    ## array 
    vcf1 <- readVcf(fl, "hg19", param=ScanVcfParam(geno="HQ")) 
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$HQ, geno(vcf2)$HQ)

    param=ScanVcfParam(geno="HQ", sample="NA00003")
    vcf1 <- readVcf(fl, "hg19", param=param) 
    suppressWarnings(writeVcf(vcf1, dest))
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$HQ, geno(vcf2)$HQ)

    ## matrix and array 
    param=ScanVcfParam(geno=c("GT", "HQ"), samples="NA00002")
    vcf1 <- readVcf(fl, "hg19", param=param)
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT)
    checkIdentical(geno(vcf1)$HQ, geno(vcf2)$HQ)
}


