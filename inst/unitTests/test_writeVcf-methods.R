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

test_alt_description_quoting = function()   # https://github.com/Bioconductor/VariantAnnotation/issues/52
{
fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(fl, genome="hg19")
tmp <- tempfile()
writeVcf(vcf, filename=tmp)
#lines = readLines(tmp) # missing from chunk above
#lines[grepl("ALT=", lines)]
require("S4Vectors")
good = new("DFrame", rownames = c("DEL", "DEL:ME:ALU", "DEL:ME:L1", 
"DUP", "DUP:TANDEM", "INS", "INS:ME:ALU", "INS:ME:L1", "INV", 
"CNV"), nrows = 10L, elementType = "ANY", elementMetadata = NULL, 
    metadata = list(), listData = list(Description = c("Deletion", 
    "Deletion of ALU element", "Deletion of L1 element", "Duplication", 
    "Tandem Duplication", "Insertion of novel sequence", "Insertion of ALU element", 
    "Insertion of L1 element", "Inversion", "Copy number variable region"
    )))
chkr = readVcf(tmp)
chk = fixed(header(readVcf(tmp)))$ALT
checkIdentical(chk, good)
}


## ---------------------------------------------------------------
## Test for GitHub issue #78: faithful round-trip via writeVcf + readVcf
## ---------------------------------------------------------------
test_writeVcf_roundtrip_issue78 <- function()
{
    ## structural.vcf has all-NA seqinfo and an existing fileDate
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    out <- tempfile(fileext=".vcf")
    on.exit(unlink(out))

    first <- readVcf(fl)
    writeVcf(first, out)
    roundtrip <- readVcf(out)

    ## Should be identical (no spurious contig lines, fileDate preserved)
    checkTrue(isTRUE(all.equal(roundtrip, first)),
        msg="structural.vcf round-trip should be perfect")

    ## ex2.vcf has real seqinfo — contig lines should still be written
    fl2 <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    out2 <- tempfile(fileext=".vcf")
    on.exit(unlink(out2), add=TRUE)

    first2 <- readVcf(fl2)
    writeVcf(first2, out2)
    roundtrip2 <- readVcf(out2)
    checkTrue(isTRUE(all.equal(roundtrip2, first2)),
        msg="ex2.vcf round-trip should be perfect")

    ## Verify fileDate is preserved (not replaced with today's date)
    lines <- readLines(out2)
    fd_line <- grep("^##fileDate=", lines, value=TRUE)
    checkIdentical(fd_line, "##fileDate=20090805")
}
