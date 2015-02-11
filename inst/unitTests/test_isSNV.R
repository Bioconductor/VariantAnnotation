vcf0 <- VCF(rowData=GRanges("chr1", IRanges(1:9, width=c(rep(1, 6), 2, 2, 3))), 
            fixed=DataFrame(
              REF=DNAStringSet(c("A", "G", "C", "T", "T", "G", "GG", "TCT", "AC")),
              ALT=DNAStringSetList("G", "A", "T", "C", c("C", "TT"), "GG", 
                                   "G", "GCG", "ACC")))
str <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(str, "")

test_isSNV_CollapsedVCF <- function() {
    checkException(isSNV(vcf), silent=TRUE)
    checkException(isInsertion(vcf), silent=TRUE)
    checkException(isDeletion(vcf), silent=TRUE)
    checkException(isIndel(vcf), silent=TRUE)
    checkException(isDelins(vcf), silent=TRUE)
    checkException(isTransition(vcf), silent=TRUE)
    checkException(isSubstitution(vcf), silent=TRUE)

    res1 <- isSNV(vcf0, singleAltOnly=TRUE)
    checkIdentical(res1, c(rep(TRUE, 4), rep(FALSE, 5)))
    res2 <- isSNV(vcf0, singleAltOnly=FALSE)
    checkIdentical(res2, c(rep(TRUE, 5), rep(FALSE, 4)))

    res1 <- isInsertion(vcf0)
    checkIdentical(res1, c(rep(FALSE, 5), TRUE, FALSE, FALSE, FALSE))
    res2 <- isDeletion(vcf0)
    checkIdentical(res2, c(rep(FALSE, 6), TRUE, FALSE, FALSE))
    res3 <- isIndel(vcf0)
    checkIdentical(res3, res1 | res2)
    res4 <- isDelins(vcf0)
    checkIdentical(res4, c(rep(FALSE, 8), TRUE))

    res1 <- isSubstitution(vcf0, singleAltOnly=TRUE)
    checkIdentical(res1, c(rep(TRUE, 4), rep(FALSE, 3), TRUE, FALSE))
    res2 <- isSubstitution(vcf0, singleAltOnly=FALSE)
    checkIdentical(res2, c(rep(TRUE, 5), rep(FALSE, 2), TRUE, FALSE))

    res1 <- isTransition(vcf0, singleAltOnly=TRUE)
    checkIdentical(res1, c(rep(TRUE, 4), rep(FALSE, 5)))
    res2 <- isTransition(vcf0, singleAltOnly=FALSE)
    checkIdentical(res2, c(rep(TRUE, 5), rep(FALSE, 4)))
}

expand <- VariantAnnotation::expand 
test_isSNV_ExpandedVCF <- function() {
    evcf <- expand(vcf)
    checkException(isSNV(evcf), silent=TRUE)
    checkException(isInsertion(evcf), silent=TRUE)
    checkException(isDeletion(evcf), silent=TRUE)
    checkException(isIndel(evcf), silent=TRUE)
    checkException(isDelins(evcf), silent=TRUE)
    checkException(isTransition(evcf), silent=TRUE)
    checkException(isSubstitution(evcf), silent=TRUE)

    evcf0 <- expand(vcf0)
    res <- isSNV(evcf0)
    checkIdentical(sum(res), 5L)

    res1 <- isInsertion(evcf0)
    checkIdentical(sum(res1), 2L)
    res2 <- isDeletion(evcf0)
    checkIdentical(sum(res2), 1L)
    res3 <- isIndel(evcf0)
    checkIdentical(sum(res3), sum(res1 + res2))
    res4 <- isDelins(evcf0)
    checkIdentical(sum(res4), 1L)

    res <- isSubstitution(evcf0)
    checkIdentical(sum(res), 6L)

    res <- isTransition(evcf0)
    checkIdentical(sum(res), 5L)
}

test_isSNV_VRanges <- function() {
    vr <- as(vcf0, "VRanges")

    res <- isSNV(vr)
    checkIdentical(sum(res), 5L)

    res1 <- isInsertion(vr)
    checkIdentical(sum(res1), 2L)
    res2 <- isDeletion(vr)
    checkIdentical(sum(res2), 1L)
    res3 <- isIndel(vr)
    checkIdentical(sum(res3), sum(res1 + res2))
    res4 <- isDelins(vr)
    checkIdentical(sum(res4), 1L)

    res <- isSubstitution(vr)
    checkIdentical(sum(res), 6L)

    res <- isTransition(vr)
    checkIdentical(sum(res), 5L)
}

test_isSNV_gvcf_format <- function() 
{
    ## ignore <NON_REF>
    fl <- system.file("unitTests", "cases", "banded_gvcf.vcf",
                      package="VariantAnnotation")

    vcf <- suppressWarnings(readVcf(fl, ""))
    checkIdentical(isSNV(vcf), rep(TRUE, nrow(vcf)))
    vr <- as(vcf, "VRanges")
    checkIdentical(isSNV(expand(vcf)), rep(TRUE, length(vr)))
}
