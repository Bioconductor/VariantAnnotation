test_snpSummary_empty <- function()
{
    cnames <- c("g00", "g01", "g11", "a0Freq", "a1Freq", 
                "HWEzscore", "HWEpvalue")
    ## no genotype
    fl <- system.file("unitTests", "cases", 
              "FORMAT_header_no_SAMPLEs.vcf", 
              package="VariantAnnotation")
    res <- suppressWarnings(snpSummary(readVcf(fl, "hg19")))
    checkTrue(all(dim(res) == c(0L, 7L)))
    checkEquals(names(res), cnames)

    ## structural
    fl <- system.file("extdata", "structural.vcf", 
              package="VariantAnnotation")
    res <- suppressWarnings(snpSummary(readVcf(fl, "hg19")))
    checkTrue(all(dim(res) == c(0L, 7L)))
    checkEquals(names(res), cnames)
}

test_snpSummary_output <- function()
{
    fl <- system.file("extdata", "ex2.vcf", 
              package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    res <- snpSummary(vcf)

    checkTrue(is.integer(res$g00))
    checkTrue(is.integer(res$g01))
    checkTrue(is.integer(res$g11))

    checkEquals(nrow(res), nrow(vcf))
    checkEquals(rownames(res), rownames(vcf))

    checkEquals(as.integer(res[1,c("g00", "g01", "g11")]), c(1L, 1L, 1L))
    checkEquals(as.integer(res[2,c("g00", "g01", "g11")]), c(2L, 1L, 0L))
    checkEquals(round(res[["a0Freq"]], 2), c(0.50, 0.83, NA, NA, NA))
    checkEquals(round(res[["a1Freq"]], 2), c(0.50, 0.17, NA, NA, NA))
}

test_snpSummary_HWEzscore <- function()
{
    ## ex2.vcf: row 1 has (1, 1, 1) -> E[g01] = 1.5 > 1 (deficit, z < 0)
    ##          row 2 has (2, 1, 0) -> E[g01] = 5/6 < 1 (excess, z > 0)
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    res <- snpSummary(readVcf(fl, "hg19"))
    checkEqualsNumeric(res$HWEzscore[1], -sqrt(1 / 3))
    checkEqualsNumeric(res$HWEzscore[2], sqrt(0.12))

    ## Synthetic VCF covering heterozygote excess (including N = 2),
    ## heterozygote deficit, and exact HWE
    gt <- rbind(
        c("0/1", "0/1", "./.", "./.", "./.", "./.", "./.", "./.", "./.", "./."),
        rep("0/1", 10L),
        c(rep("0/0", 5L), rep("1/1", 5L)),
        c("0/0", "0/1", "0/1", "1/1", rep("./.", 6L))
    )
    rownames(gt) <- paste0("snp", 1:4)
    colnames(gt) <- paste0("s", 1:10)
    rr <- GRanges("1", IRanges(1:4, width=1L),
                  paramRangeID=factor(rep(NA, 4L)))
    names(rr) <- rownames(gt)
    fixed <- DataFrame(
        REF=DNAStringSet(rep("A", 4L)),
        ALT=DNAStringSetList("G", "G", "G", "G"),
        QUAL=rep(NA_real_, 4L),
        FILTER=rep("PASS", 4L)
    )
    vcf <- VCF(rowRanges=rr, fixed=fixed, geno=SimpleList(GT=gt),
               colData=DataFrame(Samples=seq_len(ncol(gt)),
                                 row.names=colnames(gt)),
               collapsed=TRUE)
    s <- snpSummary(vcf)

    checkEqualsNumeric(s$HWEzscore, c(sqrt(2), sqrt(10), -sqrt(10), 0))
    checkEqualsNumeric(s$HWEpvalue,
                       pchisq(c(2, 10, 10, 0), df=1, lower.tail=FALSE))
}

