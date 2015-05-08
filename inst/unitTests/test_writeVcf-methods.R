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
    hd1 <- metadata(vcf1)$header
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    hd2 <- metadata(vcf2)$header
    checkTrue(names(meta(hd1)) %in% names(meta(hd2))) 
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

test_writeVcf_geno <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    dest <- tempfile()

    ## empty
    vcf1 <- readVcf(fl, "hg19", param=ScanVcfParam(geno=NA))
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkTrue(length(geno(vcf2)$GT) == 0L)

    ## Rle 
    vcf1 <- readVcf(fl, "hg19", param=ScanVcfParam(geno="GT"))
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT)
    checkIdentical(geno(header(vcf1)), geno(header(vcf2)))

    ## matrix 
    vcf1 <- readVcf(fl, "hg19", param=ScanVcfParam(geno="GT"))
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT)
    checkIdentical(geno(header(vcf1)), geno(header(vcf2)))

    param=ScanVcfParam(geno="GT", samples="NA00003")
    vcf1 <- readVcf(fl, "hg19", param=param)
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT)
    checkIdentical(geno(header(vcf1)), geno(header(vcf2)))
    checkTrue(samples(header(vcf2)) == "NA00003")

    ## array 
    vcf1 <- readVcf(fl, "hg19", param=ScanVcfParam(geno="HQ")) 
    writeVcf(vcf1, dest)
    vcf2 <-                       # FORMAT descriptors for GENO fields
        tryCatch(readVcf(dest, "hg19"), error=conditionMessage,
                 warning=conditionMessage)
    checkIdentical(geno(vcf1)$HQ, geno(vcf2)$HQ)
    checkIdentical(geno(header(vcf1)), geno(header(vcf2)))

    param=ScanVcfParam(geno="HQ", sample="NA00003")
    vcf1 <- readVcf(fl, "hg19", param=param) 
    suppressWarnings(writeVcf(vcf1, dest))
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$HQ, geno(vcf2)$HQ)
    checkIdentical(geno(header(vcf1)), geno(header(vcf2)))
    checkTrue(samples(header(vcf2)) == "NA00003")

    ## matrix and array 
    param=ScanVcfParam(geno=c("GT", "HQ"), samples="NA00002")
    vcf1 <- readVcf(fl, "hg19", param=param)
    writeVcf(vcf1, dest)
    vcf2 <- readVcf(dest, "hg19")
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT)
    checkIdentical(geno(vcf1)$HQ, geno(vcf2)$HQ)
    checkTrue(samples(header(vcf2)) == "NA00002")

    vcf3 <- vcf1
    geno(vcf3) <- geno(vcf3)[c("HQ", "GT")]
    writeVcf(vcf3, dest)
    vcf2 <- readVcf(dest, "hg19", param=param)
    checkIdentical(geno(vcf1), geno(vcf2))

    ## list
    fl <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation")
    hdr <- scanVcfHeader(fl)
    param <- ScanVcfParam(samples=samples(hdr)[1:2])
    vcf1 <- readVcf(fl, "", param=param)
    writeVcf(vcf1, dest) 
    vcf2 <- readVcf(dest, "") 
    checkIdentical(geno(vcf1)$GL, geno(vcf2)$GL)
}
