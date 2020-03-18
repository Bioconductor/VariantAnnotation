fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
scn <- scanVcf(fl)

test_FixedTypes <- function()
{
    .vcf_map_fixed <- VariantAnnotation:::.vcf_map_fixed
    exp <- exp0 <- list(rowRanges=NULL, REF=NULL,
                        ALT=list("A", character()),
                        QUAL=list("1", numeric()),
                        FILTER=list("1", character()))
    named <- names(exp)[-(1:2)]
    checkIdentical(exp, .vcf_map_fixed(character(), FALSE))
    checkIdentical(exp[1:2], .vcf_map_fixed(NA, FALSE))
    exp <- exp0
    exp[1] <- list(NULL)
    checkIdentical(exp, .vcf_map_fixed(named, FALSE))
    warn <- FALSE
    obs <- withCallingHandlers({
        .vcf_map_fixed("FOO", FALSE)
    }, warning=function(w) {
        warn <<- TRUE
        invokeRestart("muffleWarning")
    })
    checkTrue(warn)
    checkIdentical(exp[1:2], obs)
}

test_InfoTypes <- function()
{
    fmt <- info(scanVcfHeader(fl))
    info <- scn[[1]]$INFO

    checkIdentical(as.integer(c(3, 3, 2, 3, 3)), info$NS)
    checkIdentical(as.integer(c(14, 11, 10, 13, 9)), info$DP)
    checkTrue(is(info$AF,"matrix"))
    checkIdentical(c(TRUE, FALSE, TRUE, FALSE, FALSE), info$DB)
    checkIdentical(c(TRUE, rep(FALSE, 4)), info$H2)
}

test_GenoTypes <- function()
{
    fmt <- geno(scanVcfHeader(fl))
    geno <- scn[[1]]$GENO

    checkEquals(typeof(unlist(geno$GT)), "character")
    checkIdentical(lapply(geno, function(x) class(x)[1]), list(GT="matrix", GQ="matrix",
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
    checkIdentical(rep(".", 5L), info)
}

test_scanVcf_negative_Number <- function()
{
    fl <- system.file(package="VariantAnnotation", "unitTests",
                      "cases", "negative_FORMAT_Number.vcf")
    pl <- scanVcf(fl)[[1]]$GENO$PL
    checkIdentical(3L, unique(sapply(pl, length)))
    checkIdentical(761L, sum(sapply(pl, sum)))
}

test_scanVcf_crlf <- function()
{
    writeLines(readLines(fl), xx <- tempfile(), sep="\r\n")
    checkIdentical(scanVcfHeader(fl), scanVcfHeader(xx))
}

test_scanVcf_connection_n <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")

    res <- scanVcf(file(fl))
    checkIdentical(5L, length(res[[1]]$REF))

    res <- scanVcf(file(fl), n = 1000)
    checkIdentical(5L, length(res[[1]]$REF))

    res <- scanVcf(file(fl), n = 1)
    checkIdentical(1L, length(res[[1]]$REF))

    vcf <- file(fl)
    open(vcf)
    for (i in 1:10) {
        res <- scanVcf(vcf, n = 1)
        if (!length(res[[1]]$REF))
            break
    }
    close(vcf)
    checkIdentical(6L, i)
}

test_scanVcfHeader <- function()
{
    checkIdentical(VCFHeader(), scanVcfHeader())
    checkIdentical(VCFHeader(), scanVcfHeader(character()))
}

test_scanVcfHeader_VarScan <- function()
{
    fl <- system.file("unitTests", "cases", "VarScan_header.vcf",
                      package="VariantAnnotation")
    hd <- scanVcfHeader(fl)
    checkIdentical(dim(info(hd)), c(7L, 3L))
    checkIdentical(names(info(hd)), c("Number", "Type", "Description"))
    expected <- paste0("Somatic status of variant ",
        "(0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)")
    checkIdentical(info(hd)["SS", "Description"], expected)

    fl <- system.file("extdata", "ex2.vcf",
                      package="VariantAnnotation")
    hd <- scanVcfHeader(fl)
    checkIdentical(length(names(header(hd))), 11L)
    nms <- c("fileDate", "fileformat", "phasing", "reference", "source",
             "contig", "SAMPLE", "PEDIGREE")
    checkIdentical(names(meta(hd)), nms)
    checkIdentical(header(hd)$contig[["assembly"]], "B36")
}

test_scanVcfHeader_META <- function()
{
    fl <- system.file("unitTests", "cases", "meta_header.vcf",
                      package="VariantAnnotation")
    hd <- scanVcfHeader(fl)
    nms <- c("fileformat", "Tassel", "META", "SAMPLE")
    checkIdentical(names(meta(hd)), nms)
}

test_scan_row.names <- function()
{
    fl <- system.file("extdata", "chr7-sub.vcf.gz", package="VariantAnnotation")
    scn <- scanVcf(fl)[[1]]
    checkTrue(!is.null(names(scn$rowRanges)))
    scn <- scanVcf(fl, row.names=FALSE)
    checkTrue(is.null(names(scn$rowRanges)))
    param <- ScanVcfParam(which=GRanges("7", IRanges(55000723, 55000789)))
    scn <- scanVcf(fl, param=param, row.names=FALSE)
    checkTrue(is.null(names(scn$rowRanges)))
}
