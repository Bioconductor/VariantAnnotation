
test_VCF_construction <- function() {
    ## empty
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(VCF()))
    checkTrue(validObject(new("VCF")))
    checkTrue(validObject(VCF(geno=SimpleList(m1))))

    ## substance
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19") 
    checkTrue(validObject(VCF(rowData=rowData(vcf), colData=colData(vcf), 
        geno=geno(vcf), info=values(info(vcf)), 
        fixed=values(fixed(vcf))))) 
    checkIdentical("hg19", unique(genome(rowData(vcf)))) 
}

test_VCF_accessors <- function() {
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")

    ## ref
    checkTrue(class(values(ref(vcf))[["REF"]]) == "DNAStringSet")
    checkException(ref(vcf) <- NULL, silent=TRUE)
    checkException(ref(vcf) <- as.character(values(ref(vcf))[["REF"]]), 
        silent=TRUE)

    ## alt 
    checkTrue(class(values(alt(vcf))[["ALT"]]) == "DNAStringSetList")
    checkException(alt(vcf) <- NULL, silent=TRUE)
    v1 <- vcf
    alt(v1) <- CharacterList(list("A", "A", c("G", "T"), "", c("G", "GTCT"))) 
    checkTrue(validObject(v1))

    ## qual 
    checkTrue(class(values(qual(vcf))[["QUAL"]]) == "numeric")
    checkException(qual(vcf) <- NULL, silent=TRUE)
    checkException(qual(vcf) <- as.character(values(qual(vcf))[["QUAL"]]), 
        silent=TRUE)

    ## filt 
    checkTrue(class(values(filt(vcf))[["FILTER"]]) == "character")
    checkException(filt(vcf) <- NULL, silent=TRUE)
    checkException(filt(vcf) <- as.list(values(filt(vcf))[["FILTER"]]), 
        silent=TRUE)

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
    checkIdentical(ss1, ss1[TRUE,])
    checkIdentical(c(0L, ncol(ss1)), dim(ss1[FALSE,]))
    checkIdentical(ss1, ss1[,TRUE])
    checkIdentical(c(nrow(ss1), 0L), dim(ss1[,FALSE]))
    idx <- c(TRUE, FALSE)               # recycling
    ss2 <- ss1[idx,]
    checkIdentical(rowData(ss1)[idx,,drop=FALSE], rowData(ss2))
    ss2 <- ss1[,idx]
    checkIdentical(colData(ss1)[idx,,drop=FALSE], colData(ss2))
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
