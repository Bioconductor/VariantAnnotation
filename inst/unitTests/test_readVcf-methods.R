f1 <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
f2 <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
f3 <- system.file("extdata", "ex3.vcf", package="VariantAnnotation")
st <- system.file("extdata", "structural.vcf", package="VariantAnnotation")

test_readVcf_format <- function()
{
    ## arrays in geno
    vcf <- readVcf(f2, genome="hg19")
    checkTrue(class(geno(vcf)$HQ) == "array")

    ## duplicate header lines, missing INFO
    vcf <- suppressWarnings(readVcf(f3, genome="hg19"))
    checkTrue(ncol(values(info(vcf))) == 2L)
    checkTrue(class(values(alt(vcf))[["ALT"]]) == "DNAStringSetList")

    ## structural 
    vcf <- readVcf(st, genome="hg19")
    checkTrue(class(values(alt(vcf))[["ALT"]]) == "CompressedCharacterList")
    checkIdentical(values(qual(vcf))[["QUAL"]], c(NA, 6, 12, 23, 14, 11))
}

test_readVcf_ranges <- function()
{
    vcf <- readVcf(f2, "hg19")
    checkEquals(width(rowData(vcf)), width(values(ref(vcf))[["REF"]]))

    compressVcf <- bgzip(f2, tempfile())
    idx <- indexTabix(compressVcf, "vcf")
    tab <- TabixFile(compressVcf, idx)
    rd <- rowData(vcf)
    param <- ScanVcfParam(which=rd) 
    vcf_rd <- readVcf(tab, "hg19", param) 
    checkIdentical(values(info(vcf))[c("AA", "AF", "DB", "DP", "H2")], 
                   values(info(vcf_rd))[c("AA", "AF", "DB", "DP", "H2")])

    param <- ScanVcfParam(which=rd[c(3,5)]) 
    vcf_rd <- readVcf(tab, "hg19", param) 
    checkEquals(2L, dim(vcf_rd)[1])
}

test_readVcf_param <- function()
{
    gnms <- c("GT", "GQ", "DP", "HQ")
    inms <- c("NS", "DP", "AF", "DB")
 
    ## geno
    g <- gnms[2:3]
    param <- ScanVcfParam(geno=g)
    vcf <- readVcf(f2, "hg19", param)
    checkTrue(length(names(geno(vcf))) == length(g))
    checkTrue(all(names(geno(vcf)) %in% g))

    ## info 
    i <- inms[c(1,4)]
    param <- ScanVcfParam(info=i)
    vcf <- readVcf(f2, "hg19", param)
    checkTrue(ncol(values(info(vcf))) == length(i) + 1)
    checkTrue(all(i %in% names(values(info(vcf)))))

    ## geno, info combined
    param <- ScanVcfParam()
    vcf_a <- readVcf(f1, "hg19", param)
    vcf_b <- readVcf(f1, "hg19")
    checkIdentical(names(geno(vcf_a)), names(geno(vcf_b))) 
    checkIdentical(rowData(vcf_a), rowData(vcf_b))

    ## info, geno, ranges combined
    g <- gnms[1]
    i <- inms[2:3]
    rngs <- GRanges("20", IRanges(1110000, 1234600))
    param <- ScanVcfParam(geno=g, info=i, which=rngs)
    compressVcf <- bgzip(f2, tempfile())
    idx <- indexTabix(compressVcf, "vcf")
    tab <- TabixFile(compressVcf, idx)
    vcf <- readVcf(tab, "hg19",  param)
    checkTrue(all(i %in% names(values(info(vcf)))))
    checkTrue(all(names(geno(vcf)) %in% g))
    checkTrue(length(rowData(vcf)) == 3)

    ## no info, geno
    checkTrue(validObject(readVcf(f2, "hg19", ScanVcfParam(geno=NA))))
    checkTrue(validObject(readVcf(f2, "hg19", ScanVcfParam(info=NA))))
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

