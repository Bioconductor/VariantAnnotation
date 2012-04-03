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
    checkIdentical(c(TRUE, FALSE, TRUE, FALSE, FALSE), info$DB)
    checkIdentical(c(TRUE, rep(FALSE, 4)), info$H2)
}

test_GenoTypes <- function()
{
    fmt <- geno(scanVcfHeader(fl))
    geno <- scn[[1]]$GENO

    checkEquals(typeof(unlist(geno$GT)), "character")
    checkIdentical(lapply(geno, class), list(GT="matrix", GQ="matrix",
                   DP="matrix", HQ="array"))
    mat <- matrix(c(1, 3, 6, 7, 4, 8, 5, 0, 4, 2, 5, 3, 4, 2, 3), 
        nrow=5, dimnames=list(NULL, c("NA00001", "NA00002", "NA00003")))
    checkEquals(mat, geno$DP)
} 

test_scanVcf_no_FORMAT_column <- function()
{
    ## no FORMAT -- don't parse GENO
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "no_FORMAT_column.vcf")
    geno <- scanVcf(fl)[[1]]$GENO
    checkIdentical(setNames(list(), character()), geno)
}

test_scanVcf_FORMAT_header_no_SAMPLEs <- function()
{
    ## GENO tags, but no SAMPLE or actual samples
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "FORMAT_header_no_SAMPLEs.vcf")
    geno <- scanVcf(fl)[[1]]$GENO
    checkIdentical(c("GT", "DS", "GL"), names(geno))
    checkTrue(all(sapply(geno, nrow) == 5L))
    checkTrue(all(sapply(geno, ncol) == 0L))
}
    
test_scanVcf_no_INFO_header <- function()
{
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "no_INFO_header.vcf")
    info <- suppressWarnings(scanVcf(fl)[[1]]$INFO$INFO)
    checkIdentical(c(5L, 1L), dim(info))
    checkIdentical(".", unique(as.vector(info)))
}
