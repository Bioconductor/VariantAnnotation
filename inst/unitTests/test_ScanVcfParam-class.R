test_ScanVcfParam_which <- function()
{
    which <- IRangesList(seq1=IRanges(1000, 2000), 
                         seq2=IRanges(c(100, 1000), c(1000, 2000)))
    svp <- ScanVcfParam(which=which)
    checkIdentical(which, vcfWhich(svp))
    checkException(vcfWhich(svp) <- DataFrame(), silent=TRUE)
    checkException(vcfWhich(svp) <- SimpleList(), silent=TRUE)
}

test_ScanVcfParam_fixed <- function()
{
    fx <- c("GT", "ALT")
    svp <- ScanVcfParam(fixed=fx)
    checkIdentical(fx, vcfFixed(svp))
    checkException(vcfFixed(svp) <- DataFrame(), silent=TRUE)
    checkException(vcfFixed(svp) <- 1:5, silent=TRUE)

    fx <- NA_character_
    svp <- ScanVcfParam(fixed=fx)
    checkIdentical(fx, vcfFixed(svp))
}

test_ScanVcfParam_info <- function()
{
    info <- c("NS", "DP")
    svp <- ScanVcfParam(info=info)
    checkIdentical(info, vcfInfo(svp))
    checkException(vcfInfo(svp) <- DataFrame(), silent=TRUE)
    checkException(vcfInfo(svp) <- 1:5, silent=TRUE)

    info <- NA_character_
    svp <- ScanVcfParam(info=info)
    checkIdentical(info, vcfInfo(svp))
}

test_ScanVcfParam_geno <- function()
{
    geno <- c("GT", "GQ")
    svp <- ScanVcfParam(geno=geno)
    checkIdentical(geno, vcfGeno(svp))
    checkException(vcfGeno(svp) <- DataFrame(), silent=TRUE)
    checkException(vcfGeno(svp) <- 1:5, silent=TRUE)

    geno <- NA_character_
    svp <- ScanVcfParam(geno=geno)
    checkIdentical(geno, vcfGeno(svp))
}

test_ScanVcfParam_samples <- function()
{
    samp <- c("GT", "GQ")
    svp <- ScanVcfParam(samples=samp)
    checkIdentical(samp, vcfSamples(svp))
    checkException(vcfSamples(svp) <- DataFrame(), silent=TRUE)
    checkException(vcfSamples(svp) <- 1:5, silent=TRUE)
    checkException(ScanVcfParam(geno="DP", samples=NA), silent=TRUE)
    checkException(ScanVcfParam(geno=NA, samples="foo"), silent=TRUE)
}

test_ScanVcfParam_trimEmpty <- function()
{
    svp <- ScanVcfParam()
    checkIdentical(TRUE, vcfTrimEmpty(svp))
    vcfTrimEmpty(svp) <- FALSE
    checkIdentical(FALSE, vcfTrimEmpty(svp))
    checkException(vcfTrimEmpty(svp) <- "a", silent=TRUE)
    checkException(vcfTrimEmpty(svp) <- 1:5, silent=TRUE)
    ## FIXME : should this work?
    #checkException(vcfTrimEmpty(svp) <- NA, silent=TRUE)
}

test_ScanVcfParam_GRangesList <- function()
{
  fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
  compressVcf <- bgzip(fl, tempfile())
  idx <- indexTabix(compressVcf, "vcf")
  tab <- TabixFile(compressVcf, idx)
  gr1 <- GRanges("1", IRanges(13219, 2827695, name="regionA"))
  gr2 <- GRanges(rep("2", 2), 
      IRanges(c(321680, 14477080), c(321689, 14477090),
      name=c("regionB", "regionC")))
  grl <- GRangesList("1"=gr1, "2"=gr2)
  vcf_grl <- readVcf(tab, "hg19", grl)
  gr <- unlist(grl, use.names=FALSE)
  vcf_gr <- readVcf(tab, "hg19", gr)
  checkIdentical(rowRanges(vcf_grl), rowRanges(vcf_gr))
}
