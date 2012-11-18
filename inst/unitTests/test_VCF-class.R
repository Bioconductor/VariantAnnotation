
test_VCF_construction <- function() {
    ## empty
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(VCF()))
    checkTrue(validObject(new("CollapsedVCF")))
    checkTrue(validObject(new("ExpandedVCF")))
    checkTrue(validObject(VCF(geno=SimpleList(m1))))

    ## substance
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    target <- readVcf(fl, genome="hg19")
    current <- VCF(rowData=rowData(target), colData=colData(target),
                   geno=geno(target), info=values(info(target)),
                   fixed=values(fixed(target))) 
    checkTrue(validObject(current)) 
    checkIdentical("hg19", unique(genome(rowData(target)))) 
}

test_VCF_accessors <- function() {
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")

    ## ref
    checkTrue(class(ref(vcf)) == "DNAStringSet")
    checkException(ref(vcf) <- NULL, silent=TRUE)
    checkException(ref(vcf) <- as.character(ref(vcf)), 
        silent=TRUE)

    ## alt 
    checkTrue(class(alt(vcf)) == "DNAStringSetList")
    checkException(alt(vcf) <- NULL, silent=TRUE)
    v1 <- vcf
    alt(v1) <- CharacterList(list("A", "A", c("G", "T"), "", c("G", "GTCT"))) 
    checkTrue(validObject(v1))

    ## qual 
    checkTrue(class(qual(vcf)) == "numeric")
    checkException(qual(vcf) <- NULL, silent=TRUE)
    checkException(qual(vcf) <- as.character(qual(vcf)), silent=TRUE)

    ## filt 
    checkTrue(class(filt(vcf)) == "character")
    checkException(filt(vcf) <- NULL, silent=TRUE)
    checkException(filt(vcf) <- as.list(filt(vcf)), silent=TRUE)

    ## fixed 
    checkTrue(class(values(fixed(vcf))) == "DataFrame")
    checkException(fixed(vcf) <- NULL, silent=TRUE)
    checkException(fixed(vcf) <- as.matrix(values(fixed(vcf))), 
        silent=TRUE)

    DF <- DataFrame(REF=DNAStringSet(c("A", "C", "T", "T", "C")),
                    ALT=CharacterList(list("C", "G", "A", "G", "G")),
                    QUAL=c(10, 20, NA, 10, 10),
                    FILTER=c("pass", "pass", "q10", "q10", "pass"))
    df <- DF
    names(df) <- c("reference", "ALT", "QUAL", "FILTER")
    checkException(fixed(vcf) <- df, silent=TRUE)
    df <- DF
    names(df) <- c("REF", "ALT", "qual", "FILTER")
    checkException(fixed(vcf) <- df, silent=TRUE)
    df <- DF
    df$QUAL <- as.character(df$QUAL)
    checkException(fixed(vcf) <- df, silent=TRUE)
    df <- DF
    df$ALT <- unlist(df$ALT, use.names=FALSE)
    checkException(fixed(vcf) <- df, silent=TRUE)
    df <- DF
    df$REF <- as.character(df$REF)
    checkException(fixed(vcf) <- df, silent=TRUE)

    ## info
    checkTrue(class(values(info(vcf))) == "DataFrame")
    checkException(info(vcf) <- NULL, silent=TRUE)
    v1 <- vcf 
    info(v1) <- DataFrame()

    checkTrue(class(values(info(vcf))[["AF"]]) == "CompressedNumericList")
    AF <- NumericList(0.5, 0.017, c(0.333,0.667), NA, NA)
    checkIdentical(values(info(vcf))[["AF"]], AF) 

    ## geno
    checkTrue(class(geno(vcf)) == "SimpleList")
    checkException(geno(vcf) <- NULL, silent=TRUE)
    checkIdentical(geno(vcf)$CNQ, geno(vcf)[["CNQ"]])
    checkIdentical(geno(vcf)[["CN"]], assays(vcf)[["CN"]]) 
}

test_VCF_subset <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")

    ## numeric
    ss1 <- vcf[2:3,]
    checkIdentical(c(2L, ncol(vcf)), dim(ss1))
    ss1 <- vcf[,2]
    checkIdentical(c(nrow(vcf), 1L), dim(ss1))
    checkIdentical(rownames(ss1), names(info(ss1)))
    ss1 <- vcf[2:3, 2]
    checkIdentical(c(2L, 1L), dim(ss1))
    checkIdentical(rownames(ss1), names(info(ss1)))

    ## character
    ss1 <- vcf 
    dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))], letters[seq_len(ncol(ss1))]) 
    ridx <- c("B", "C") 
    checkIdentical(rowData(ss1[ridx,]), rowData(ss1)[ridx,]) 
    checkIdentical(rowData(ss1["C",]), rowData(ss1)["C",,drop=FALSE])
    checkException(ss1[LETTERS,], "i-index out of bounds", TRUE)
    cidx <- "b"
    checkIdentical(colData(ss1[,cidx]), colData(ss1)[cidx,,drop=FALSE])
    checkException(ss1[,letters], "j-index out of bounds", TRUE)

    ## logical
    ss1 <- vcf 
    dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                          letters[seq_len(ncol(ss1))])
    checkIdentical(fixed(ss1), fixed(ss1[TRUE,]))
    checkIdentical(info(ss1), info(ss1[TRUE,]))
    checkIdentical(c(0L, ncol(ss1)), dim(ss1[FALSE,]))
    checkIdentical(fixed(ss1), fixed(ss1[,TRUE]))
    checkIdentical(info(ss1), info(ss1[,TRUE]))
    checkIdentical(c(nrow(ss1), 0L), dim(ss1[,FALSE]))
    idx <- c(TRUE, FALSE)               # recycling
    ss2 <- ss1[idx,]
    checkIdentical(rowData(ss1)[idx,,drop=FALSE], rowData(ss2))
    ss2 <- ss1[,idx]
    checkIdentical(colData(ss1)[idx,,drop=FALSE], colData(ss2))

    ## 0 columns
    vcf <- VCF(rowData=GRanges("chr1", IRanges(1:10, width=1)))
    checkIdentical(dim(vcf[1:5, ]), c(5L, 0L))
    ## 0 rows 
    vcf <- VCF(colData=DataFrame(samples=1:10))
    checkIdentical(dim(vcf[ ,1:5]), c(0L, 5L))
}

test_VCF_subset_empty_slots <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19", param=ScanVcfParam(info=NA, fixed=NA))
    checkTrue(names(mcols(vcf)) == "paramRangeID") 
    checkTrue(nrow(vcf[2:4]) == 3L) 
    checkTrue(ncol(vcf[,1:2]) == 2L) 

    vcf <- readVcf(fl, "hg19", param=ScanVcfParam(geno=NA))
    checkTrue(length(names(geno(vcf))) == 0L) 
    checkTrue(ncol(vcf) == 0L) 
    checkException(vcf[,1:5], silent=TRUE) 
    checkTrue(nrow(vcf[2:4]) == 3L)
}

test_VCF_subsetassign <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")
    ss1 <- vcf 
    ss1[1:2,] <- ss1[2:1,]
    checkIdentical(rowData(vcf)[2:1,], rowData(ss1)[1:2,])
    checkIdentical(rowData(vcf[-(1:2),]), rowData(ss1)[-(1:2),])
    checkIdentical(colData(vcf), colData(ss1))
    checkIdentical(c(exptData(vcf), exptData(vcf)), exptData(ss1))

    ss1 <- vcf 
    ss1[,1:2] <- ss1[,2:1,drop=FALSE]
    checkIdentical(colData(vcf)[2:1,,drop=FALSE],
                   colData(ss1)[1:2,,drop=FALSE])
    checkIdentical(colData(vcf)[-(1:2),,drop=FALSE],
                   colData(ss1)[-(1:2),,drop=FALSE])
    checkIdentical(rowData(vcf), rowData(ss1))
    checkIdentical(c(exptData(vcf), exptData(vcf)), exptData(ss1))
}

test_VCF_seqlevels <- function()
{
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")
    seqlev <- seqlevels(vcf)
    checkIdentical(seqlev, c("1", "2", "3", "4"))

    vcf2 <- renameSeqlevels(vcf, c("1"="chr1"))
    checkIdentical(seqlevels(vcf2), c("chr1", "2", "3", "4"))
 
    vcf3 <- keepSeqlevels(vcf, "3")
    checkIdentical(seqlevels(vcf3), "3")
}
