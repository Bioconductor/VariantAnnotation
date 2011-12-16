f1 <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
f2 <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
f3 <- system.file("extdata", "structural.vcf", package="VariantAnnotation")

test_readVcf_ranges <- function()
{
    vcf <- readVcf(f1, "hg19")
    checkEquals(width(rowData(vcf)), width(values(rowData(vcf))[["REF"]]))
    checkEquals(scanVcf(f1)[[1]]$POS, start(rowData(vcf))) 
}

test_readVcf_metadata <- function()
{
    vcf <- readVcf(f1, "hg19")
    checkTrue(all(rownames(exptData(vcf)[["HEADER"]][["FORMAT"]]) %in% 
              names(assays(vcf)))) 
    checkIdentical("hg19", unique(genome(rowData(vcf)))) 
}

test_readVcf_accessors <- function()
{
    vcf <- readVcf(f1, "hg19")
    checkTrue(class(vcf) == "SummarizedExperiment")
    checkTrue(class(rowData(vcf)) == "GRanges")
    checkTrue(class(assays(vcf)) == "SimpleList")
    checkTrue(class(exptData(vcf)) == "SimpleList")
    checkTrue(class(colData(vcf)) == "DataFrame")
}

test_readVcf_formats <- function()
{
    ## arrays in @assays 
    vcf <- readVcf(f2, "hg19")
    checkTrue(class(assays(vcf)$HQ) == "array")

    ## structural 
    vcf <- readVcf(f3, "hg19")
    checkTrue(class(rowData(vcf)) == "GRanges")
}

test_readVcf_param <- function()
{
    ## geno
    param <- ScanVcfParam(geno=c("DP", "GQ"))
    vcf <- readVcf(f2, "hg19", param=param)
    checkTrue(length(names(assays(vcf))) == 2)
    checkTrue(all(names(assays(vcf)) %in% c("DP", "GQ")))

    param <- ScanVcfParam(geno=c("DP", "DL"))
    vcf <- readVcf(f2, "hg19", param=param)
    checkTrue(names(assays(vcf)) == "DP")

    ## info 
    param <- ScanVcfParam(info="NS")
    vcf <- readVcf(f2, "hg19", param=param)
    #checkTrue("NS.1" %in% colnames(values(rowData(vcf))))

    ## geno, info
    param <- ScanVcfParam()
    vcf_a <- readVcf(f1, "hg19", param=param)
    vcf_b <- readVcf(f1, "hg19")
    checkIdentical(names(assays(vcf_a)), names(assays(vcf_b))) 
    checkIdentical(rowData(vcf_a), rowData(vcf_b))

    ## info, geno, ranges
    rngs <- GRanges("20", IRanges(1110000, 1234600))
    param <- ScanVcfParam(geno="HQ", info="AF", which=rngs)
    compressVcf <- bgzip(f2, tempfile())
    idx <- indexTabix(compressVcf, "vcf")
    tab <- TabixFile(compressVcf, idx)
    vcf <- readVcf(tab, "hg19", param=param)
    checkTrue("AF" %in% colnames(values(rowData(vcf))))
    checkTrue(all(names(assays(vcf)) %in% "HQ"))
    checkTrue(length(rowData(vcf)) == 3)
}


test_readVcf_tabix <- function()
{
    ## variant at position 101558 
    param1 <- GRanges(seqnames="16", ranges=IRanges(start=103466, end=103476))
    param2 <- GRanges(seqnames="16", ranges=IRanges(start=103476, end=103476))
    param3 <- GRanges(seqnames="16", ranges=IRanges(start=103476, end=103486))
    cmp <- bgzip(f1, tempfile())
    idx <- indexTabix(cmp, "vcf")
    tbx <- TabixFile(cmp, idx)
 
    scn1 <- scanTabix(tbx, param=param1)
    scn2 <- scanTabix(tbx, param=param2)
    scn3 <- scanTabix(tbx, param=param3)
    #checkIdentical(scn1, scn2)
    #checkIdentical(scn2, scn3)
    #checkIdentical(scn1, scn3)
}

