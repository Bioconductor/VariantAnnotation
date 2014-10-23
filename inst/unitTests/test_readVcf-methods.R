test_readVcf_format <- function()
{
    ## arrays in geno
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    checkTrue(class(geno(vcf)$HQ) == "array")

    ## missing QUAL, FILTER, INFO 
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "no_INFO_header.vcf")
    vcf <- suppressWarnings(readVcf(fl, "hg19"))
    checkTrue("DNAStringSetList" == class(alt(vcf)))
    checkTrue("numeric" == class(qual(vcf)))
    checkTrue("character" == class(filt(vcf)))

    ## structural 
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    checkTrue(class(alt(vcf)) == "CompressedCharacterList")
    checkIdentical(qual(vcf), c(6, NA, 6, 12, 23, 14, 11))
}

test_readVcf_unspecified_INFO_FORMAT <- function()
{
    ## As of 1.7.32, warnings are no longer thrown for 'INFO'
    ## and 'FORMAT' fields in the data but with no header entry.
    ## Fields are silently skipped.
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "unspecified_INFO_FORMAT_fields.vcf")
    vcf <- readVcf(fl, "hg19")
    ## columns immediately after XX entries
    exp <- c(14L, 11L, 10L, 13L, 9L)
    checkIdentical(exp, info(vcf)$DP)
    exp <- NumericList(0.5, 0.017, c(0.333, 0.667), NA, c(NA, NA))
    checkIdentical(exp, info(vcf)$AF)

    ## columns immediately after FORMAT entries
    exp <- c("0|0", "0/1", "0|0", "0/2", "0/0", "1/1")
    checkIdentical(exp, as.vector(geno(vcf)$GT[4:5,]))
    exp <- c(56L, NA, 51L, NA, NA, NA, 60L, NA, 51L, NA, NA, NA)
    checkIdentical(exp, as.vector(geno(vcf)$HQ[4:5,,]))
}

test_readVcf_fewer_FORMAT_than_GENO <- function() {
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "fewer-FORMAT-than-GENO.vcf")
    exp <- "record 2 sample H1993: fewer FORMAT fields than GENO fields"
    obs <- tryCatch(readVcf(fl, "hg19"), warning=conditionMessage)
    checkIdentical(exp, obs)
}

test_readVcf_missing_FORMAT_metadata_elt <- function() {
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "missing-FORMAT-metadata-elt.vcf")
    vcf <- readVcf(fl, "hg19")
    exp <- c(1L, NA_integer_, 1L) 
    obs <- unlist(geno(vcf)$AP, use.names=FALSE)
    checkIdentical(exp,  obs)
}

test_readVcf_no_GENO_row <- function() {
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "no_GENO_row.vcf")
    exp <- c(582L, NA, 584L, 583L, NA, 585L)
    vcf <- readVcf(fl, "hg19")
    checkIdentical(exp, as.vector(geno(vcf)[["DP"]]))
}

test_readVcf_ranges <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    checkEquals(width(rowData(vcf)), width(ref(vcf)))

    compressVcf <- bgzip(fl, tempfile())
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
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    gnms <- c("GT", "GQ", "DP", "HQ")
    inms <- c("NS", "DP", "AF", "DB")
    snms <- c("NA00001", "NA00002", "NA00003")

    ## samples
    samp <- snms[3:2]
    vcf <- readVcf(fl, "hg19", param=ScanVcfParam(samples=samp))
    checkTrue(ncol(vcf) == 2L)
    checkIdentical(colnames(vcf), snms[3:2])

    param <- ScanVcfParam(geno=c("GT", "HQ"), samples=snms)
    vcf1 <- readVcf(fl, "hg19")
    vcf2 <- readVcf(fl, "hg19", param=param)
    checkIdentical(geno(vcf1)$GT, geno(vcf2)$GT) 
    checkIdentical(geno(vcf1)$HQ, geno(vcf2)$HQ)
 
    samp <- snms[3] 
    param <- ScanVcfParam(geno=c("GT", "HQ"), samples=samp)
    vcf3 <- readVcf(fl, "hg19", param=param)
    current <- c(geno(vcf3)$GT)
    checkIdentical(unname(geno(vcf1)$GT[,samp]), current)
    current <- c(geno(vcf3)$HQ)
    checkIdentical(c(geno(vcf1)$HQ[,samp,]), current)

    ## geno
    g <- gnms[3:2]
    param <- ScanVcfParam(geno=g)
    vcf <- readVcf(fl, "hg19", param)
    checkTrue(length(names(geno(vcf))) == length(g))
    checkIdentical(names(geno(vcf)), gnms[3:2])

    fl <- system.file("extdata", "chr22.vcf.gz",
                      package="VariantAnnotation")
    param1 <- ScanVcfParam(which=GRanges("22", IRanges(5e7, 50302629)))
    vcf1 <- readVcf(fl, "hg19", param=param1)
    param2 <- ScanVcfParam(which=rowData(vcf1)[1:10])
    vcf2 <- readVcf(fl, "hg19", param=param2)
    checkIdentical(geno(vcf1)$GT[1:10, ], geno(vcf2)$GT)

    ## info 
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    i <- inms[c(4,1)]
    param <- ScanVcfParam(info=i)
    vcf <- readVcf(fl, "hg19", param)
    checkTrue(ncol(info(vcf)) == length(i))
    checkIdentical(colnames(info(vcf)), inms[c(4,1)])

    ## geno, info
    param <- ScanVcfParam()
    vcf1 <- readVcf(fl, "hg19", param)
    vcf2 <- readVcf(fl, "hg19")
    checkIdentical(names(geno(vcf1)), names(geno(vcf2))) 
    checkIdentical(rowData(vcf1), rowData(vcf2))

    ## info, geno, ranges, samples
    g <- gnms[1]
    i <- inms[2:3]
    s <- snms[3]
    rngs <- GRanges("20", IRanges(1110000, 1234600))
    param <- ScanVcfParam(geno=g, info=i, which=rngs, samples=s)
    compressVcf <- bgzip(fl, tempfile())
    idx <- indexTabix(compressVcf, "vcf")
    tab <- TabixFile(compressVcf, idx)
    vcf <- readVcf(tab, "hg19",  param)
    checkTrue(all(i %in% colnames(info(vcf))))
    checkTrue(all(names(geno(vcf)) %in% g))
    checkTrue(length(rowData(vcf)) == 3)
    checkTrue(ncol(vcf) == 1L)
    checkTrue(colnames(vcf) == s)

    ## no info, geno, samples
    checkTrue(validObject(readVcf(fl, "hg19", ScanVcfParam(geno=NA))))
    checkTrue(validObject(readVcf(fl, "hg19", ScanVcfParam(info=NA))))
    obs <-                              # no warnings on samples=NA
        tryCatch(readVcf(fl, "hg19", ScanVcfParam(samples=NA)),
                 warning=conditionMessage)
    checkTrue(is(obs, "VCF") && validObject(obs))
}

test_readVcf_genome <- function()
{
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    checkException(readVcf(fl, genome=2), silent=TRUE)

    vcf0 <- readVcf(fl, "test")
    si0 <- seqinfo(vcf0)

    ## no param 
    si1 <- si0 
    isCircular(si1) <- rep(TRUE, 4)
    vcf1 <- readVcf(fl, si1)
    checkIdentical(si1, seqinfo(vcf1))
    si2 <- si0[c("1", "2")]
    vcf2 <- readVcf(fl, si2)
    checkIdentical(unname(genome(vcf2)), c("test", "test", NA, NA))
    si3 <- Seqinfo(as.character(1:5), NA, NA, "test")
    vcf3 <- readVcf(fl, si3)
    checkIdentical(merge(si3, si0), seqinfo(vcf3))

    ## param
    bgz <- bgzip(fl, tempfile())
    idx <- indexTabix(bgz, "vcf")
    tab <- TabixFile(bgz, idx)

    param <- GRanges(c("2", "3"), IRanges(c(321682, 12665100), width=1))
    si4 <- Seqinfo(c("2", "3"), NA, NA) 
    vcf4 <- readVcf(tab, si4, param=param) 
    checkIdentical(si4, seqinfo(rowData(vcf4)))
    si5 <- Seqinfo("2", NA, NA) 
    checkException(readVcf(tab, si5, param=param), silent=TRUE) 
    si6 <- Seqinfo(c("1", "2", "3"), NA, NA) 
    vcf6 <- readVcf(tab, si6, param=param)
    checkIdentical(merge(si6, si4), seqinfo(vcf6))
}

test_readVcf_tabix <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")

    param1 <- GRanges(seqnames="20", ranges=IRanges(start=17320, end=17330))
    param2 <- GRanges(seqnames="20", ranges=IRanges(start=17330, end=17330))
    param3 <- GRanges(seqnames="20", ranges=IRanges(start=17330, end=17340))
    cmp <- bgzip(fl, tempfile())
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

test_readGT <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    GT <- readGT(fl, nucleotides=TRUE)
    checkIdentical(colnames(GT), c("NA00001", "NA00002", "NA00003")) 
    checkIdentical(unname(GT[1,]), c("G|G", "A|G", "A/A")) 
    checkIdentical(unname(GT[2,]), c("T|T", "T|A", "T/T"))
    checkIdentical(unname(GT[5,]), c("GTC/G", "GTC/GTCT", "G/G"))

    fl <- system.file("unitTests", "cases", "no_INFO_header.vcf", 
                      package="VariantAnnotation")
    GT <- readGT(fl, nucleotides=TRUE)
    checkIdentical(unname(GT[,1]), rep(".", 5)) 
}
