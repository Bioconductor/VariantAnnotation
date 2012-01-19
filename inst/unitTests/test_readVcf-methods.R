f1 <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
f2 <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
f3 <- system.file("extdata", "ex3.vcf", package="VariantAnnotation")
st <- system.file("extdata", "structural.vcf", package="VariantAnnotation")

test_readVcf_format <- function()
{
    ## arrays in geno
    vcf <- readVcf(f2, "hg19")
    checkTrue(class(geno(vcf)$HQ) == "array")

    ## duplicate header lines, missing INFO
    vcf <- suppressWarnings(readVcf(f3, "hg19"))
    checkTrue(ncol(values(info(vcf))) == 1L)
    checkTrue(class(values(alt(vcf))[["ALT"]]) == "DNAStringSetList")

    ## structural 
    vcf <- readVcf(st, "hg19")
    checkTrue(class(values(alt(vcf))[["ALT"]]) == "CompressedCharacterList")
    checkIdentical(values(qual(vcf))[["QUAL"]], c(NA, 6, 12, 23, 14, 11))
}

test_readVcf_accessors <- function()
{
    vcf <- readVcf(f2, "hg19")
    checkTrue(class(rowData(vcf)) == "GRanges")
    checkTrue(class(exptData(vcf)) == "SimpleList")
    checkTrue(class(colData(vcf)) == "DataFrame")
    checkTrue(class(assays(vcf)) == "SimpleList")
    checkTrue(class(geno(vcf)) == "SimpleList")
    checkIdentical(assays(vcf), geno(vcf))

    checkTrue(class(values(info(vcf))) == "DataFrame")
    checkTrue(class(values(info(vcf))[["AF"]]) == "CompressedNumericList")
    AF <- NumericList(0.5, 0.017, c(0.333,0.667), NA, NA)
    checkIdentical(values(info(vcf))[["AF"]], AF) 

    checkIdentical("hg19", unique(genome(rowData(vcf)))) 
}

test_readVcf_ranges <- function()
{
    vcf <- readVcf(f1, "hg19")
    checkEquals(width(rowData(vcf)), width(values(ref(vcf))[["REF"]]))
    checkEquals(scanVcf(f1)[[1]]$POS, start(rowData(vcf))) 
}

test_readVcf_param <- function()
{
    gnms <- c("GT", "GQ", "DP", "HQ")
    inms <- c("NS", "DP", "AF", "DB")
 
    ## geno
    g <- gnms[2:3]
    param <- ScanVcfParam(geno=g)
    vcf <- readVcf(f2, "hg19", param=param)
    checkTrue(length(names(geno(vcf))) == length(g))
    checkTrue(all(names(geno(vcf)) %in% g))

    ## info 
    i <- inms[c(1,4)]
    param <- ScanVcfParam(info=i)
    vcf <- readVcf(f2, "hg19", param=param)
    checkTrue(ncol(values(info(vcf))) == length(i))
    checkTrue(all(names(values(info(vcf))) %in% i))

    ## geno, info
    param <- ScanVcfParam()
    vcf_a <- readVcf(f1, "hg19", param=param)
    vcf_b <- readVcf(f1, "hg19")
    checkIdentical(names(geno(vcf_a)), names(geno(vcf_b))) 
    checkIdentical(rowData(vcf_a), rowData(vcf_b))

    ## info, geno, ranges
    g <- gnms[1]
    i <- inms[2:3]
    rngs <- GRanges("20", IRanges(1110000, 1234600))
    param <- ScanVcfParam(geno=g, info=i, which=rngs)
    compressVcf <- bgzip(f2, tempfile())
    idx <- indexTabix(compressVcf, "vcf")
    tab <- TabixFile(compressVcf, idx)
    vcf <- readVcf(tab, "hg19", param=param)
    checkTrue(all(names(values(info(vcf))) %in% i))
    checkTrue(all(names(geno(vcf)) %in% g))
    checkTrue(length(rowData(vcf)) == 3)
}

test_readVcf_tabix <- function()
{
    param1 <- GRanges(seqnames="16", ranges=IRanges(start=103466, end=103476))
    param2 <- GRanges(seqnames="16", ranges=IRanges(start=103476, end=103476))
    param3 <- GRanges(seqnames="16", ranges=IRanges(start=103476, end=103486))
    cmp <- bgzip(f1, tempfile())
    idx <- indexTabix(cmp, "vcf")
    tbx <- TabixFile(cmp, idx)
 
    scn1 <- scanTabix(tbx, param=param1)
    scn2 <- scanTabix(tbx, param=param2)
    scn3 <- scanTabix(tbx, param=param3)
    names(scn1) <- names(scn2) <- names(scn3) <- NULL
    checkIdentical(scn1, scn2)
    checkIdentical(scn2, scn3)
    checkIdentical(scn1, scn3)
}

