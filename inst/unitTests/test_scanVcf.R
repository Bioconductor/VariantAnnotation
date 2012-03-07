fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
scn <- scanVcf(fl)


test_FixedTypes <- function()
{
    .vcf_fixed <- VariantAnnotation:::.vcf_fixed
    exp <- exp0 <- list(rowData=NULL, REF=NULL, ALT=character(),
                        QUAL=numeric(), FILTER=character())
    checkIdentical(exp, .vcf_fixed(character()))
    exp[] <- list(NULL)
    checkIdentical(exp, .vcf_fixed(NA))
    exp <- exp0
    exp[1] <- list(NULL)
    checkIdentical(exp, .vcf_fixed(names(exp)[-(1:2)]))
    warn <- FALSE
    exp[] <- list(NULL)
    obs <- withCallingHandlers({
        .vcf_fixed("FOO")
    }, warning=function(w) {
        warn <<- TRUE
        invokeRestart("muffleWarning")
    })
    checkTrue(warn)
    checkIdentical(exp, obs)
}

test_InfoTypes <- function()
{
    fmt <- info(scanVcfHeader(fl))
    info <- scn[[1]]$INFO 

    checkIdentical(as.integer(c(3, 3, 2, 3, 3)), info$NS)
    checkIdentical(as.integer(c(14, 11, 10, 13, 9)), info$DP)
    checkEquals(class(info$AF), "matrix")
    checkIdentical(c(TRUE, FALSE, FALSE, FALSE, FALSE), info$DB)
    checkIdentical(rep(FALSE, 5), info$H2)
}

test_GenoTypes <- function()
{
    fmt <- geno(scanVcfHeader(fl))
    geno <- scn[[1]]$GENO

    checkEquals(typeof(unlist(geno$GT)), "character")
    checkIdentical(lapply(geno, class), list(GT="matrix", GQ="matrix",
                   DP="matrix", HQ="array"))
    mat <- matrix(c(1, 3, 6, 7, 4, 8, 5, 0, 4, 2), nrow=5, dimnames=list(NULL,
        c("NA001", "NA002")))
    checkEquals(mat, geno$DP)
} 
