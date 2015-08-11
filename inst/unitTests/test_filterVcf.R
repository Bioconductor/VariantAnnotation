quiet <- suppressMessages

test_filterVcf_TabixFile <- function()
{
    fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
    tbx <- TabixFile(fl, yieldSize=5000)
    dest <- tempfile()
    filt <- FilterRules(list(fun=function(...) TRUE))
    ans <- quiet(filterVcf(tbx, "hg19", dest, filters=filt))
 
    checkIdentical(dest, ans)
    vcf0 <- readVcf(fl, "hg19")
    vcf1 <- readVcf(dest, "hg19")
    checkIdentical(dim(vcf0), dim(vcf1))

    ## with ranges 
    param <- ScanVcfParam(which=GRanges("22", IRanges(50301340, width=10000)))
    ans1 <- quiet(filterVcf(tbx, "", tempfile(), filters=filt, param=param))
    ans2 <- quiet(filterVcf(fl, "", tempfile(), filters=filt, param=param))
    vcf1 <- readVcf(ans1, "")
    vcf2 <- readVcf(ans2, "")
    checkIdentical(rowRanges(vcf1), rowRanges(vcf2))
}

test_filterVcf_filter <- function()
{
    fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
    tbx <- TabixFile(fl, yieldSize=5000)

    filt <- FilterRules(list(filt1=function(x) {
        rowSums(geno(x)$DS > 0.5) > 0
    }))
    dest <- tempfile()

    ans <- filterVcf(tbx, "hg19", dest, filters=filt, verbose=FALSE)
    vcf0 <- subsetByFilter(readVcf(fl, "hg19"), filt)
    vcf1 <- readVcf(ans, "hg19")
    checkIdentical(dim(vcf0), dim(vcf1))
}

test_filterVcf_prefilter_only <- function()
{
    fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
    tbx <- TabixFile(fl, yieldSize=5000)

    filt <- FilterRules(list(filt1=function(x) {
        grepl("LOWCOV", x, fixed=TRUE)
    }))
    dest <- tempfile()

    ans <- filterVcf(tbx, "hg19", dest, prefilters=filt, verbose=FALSE)

    vcf <- readVcf(fl, "hg19")
    idx <- any(info(vcf)$SNPSOURCE == "LOWCOV")
    vcf0 <- vcf[!is.na(idx) & idx,]
    vcf1 <- readVcf(ans, "hg19")
    checkIdentical(dim(vcf0), dim(vcf1))
}

test_filterOnSomaticStatus <- function()
{
    f <- system.file("extdata", "h1187-10k.vcf.gz", package="VariantAnnotation")
    vcf <- readVcf(f, "hg19")
    somaticStatus <- as.list(table(info(vcf)$SS))
    checkEquals(somaticStatus, list(Germline=103, LOH=1, Somatic=4))

    somaticStatusGermlineFilter <- function(x){
        !is.na(info(x)$SS) & info(x)$SS=="Germline"
        }
 
        # filter the in-memory data structure
    filters <- FilterRules(list(filter.1=somaticStatusGermlineFilter))
    checkEquals(dim(subsetByFilter(vcf, filters)), c(103, 2))
}

test_prefilterOnSomaticStatus <- function()
{
    f <- system.file("extdata", "h1187-10k.vcf.gz", package="VariantAnnotation")
    tabix.file <- TabixFile(f, yieldSize=1000)

    isGermline=function(x) {
        grepl("Germline", x, fixed=TRUE)
        }

    prefilteringFunctions <- FilterRules(list(isGermline=isGermline))
    filtered.filename <- filterVcf(tabix.file, "hg19", tempfile(),
                                   prefilters=prefilteringFunctions,
                                   verbose=FALSE)

    checkEquals(nrow(readVcf(filtered.filename, "hg19")), 103)
} 

test_filterOnSnps <- function()
{
    f <- system.file("extdata", "h1187-10k.vcf.gz", package="VariantAnnotation")
    tabix.file <- TabixFile(f, yieldSize=1000)

       # filter only on snp
    isSnp=function(x) {
        refSnp <- nchar(ref(x)) == 1L
        a <- alt(x)
        altSnp <- elementLengths(a) == 1L
        ai <- unlist(a[altSnp])    # all length 1, so unlisting is 1:1 map
        altSnp[altSnp] <- nchar(ai) == 1L & (ai %in% c("A", "C", "G", "T"))
        refSnp & altSnp
        }
 
    filteringFunctions <- FilterRules(list(isSnp=isSnp))
    filtered.filename <- filterVcf(tabix.file, "hg19", tempfile(),
                                   filters=filteringFunctions,
                                   verbose=FALSE)
      # now check the results
    vcf.germline <- readVcf(filtered.filename, "hg19")   # 
    checkEquals(nrow(vcf.germline), 186)
}

test_prefilterOnSomaticStatusThenFilterOnSnps <- function()
{
    f <- system.file("extdata", "h1187-10k.vcf.gz", package="VariantAnnotation")
    tabix.file <- TabixFile(f, yieldSize=1000)

    isGermline=function(x) {
        grepl("Germline", x, fixed=TRUE)
        }

    prefilteringFunctions <- FilterRules(list(isGermline=isGermline))

    isSnp=function(x) {
        refSnp <- nchar(ref(x)) == 1L
        a <- alt(x)
        altSnp <- elementLengths(a) == 1L
        ai <- unlist(a[altSnp])    # all length 1, so unlisting is 1:1 map
        altSnp[altSnp] <- nchar(ai) == 1L & (ai %in% c("A", "C", "G", "T"))
        refSnp & altSnp
        }
 
    filteringFunctions <- FilterRules(list(isSnp=isSnp))
 
        # now filter on both.  should be 98 rows
    filtered.filename <- filterVcf(tabix.file, "hg19", tempfile(),
                                   prefilters=prefilteringFunctions,
                                   filters=filteringFunctions,
                                   verbose=FALSE)

      # now check the results
    vcf.germline.snp <- readVcf(filtered.filename, "hg19")   # 
    checkEquals(nrow(vcf.germline.snp), 98)
}
