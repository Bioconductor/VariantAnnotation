
test_VCF_construction <- function() {
    ## empty
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(VCF()))
    checkTrue(validObject(new("VCF")))
    checkTrue(validObject(VCF(assays=SimpleList(m1))))

    ## substance
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    dat <- readVcf(fl, "hg19") 
    checkTrue(validObject(VCF(rowData=rowData(dat), colData=colData(dat), 
        geno=geno(dat), info=values(info(dat)), 
        fixedFields=values(fixedFields(dat))))) 
    checkException(VCF(rowData=rowData(dat), colData=colData(dat), 
        geno=geno(dat), info=values(info(dat))), silent=TRUE) 
    checkException(VCF(rowData=rowData(dat), colData=colData(dat), 
        geno=geno(dat), fixedFields=values(fixedFields(dat))), silent=TRUE) 
}

test_VCF_accessors <- function() {
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")

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

    ## fixedFields 
    checkTrue(class(values(fixedFields(vcf))) == "DataFrame")
    checkException(fixedFields(vcf) <- NULL, silent=TRUE)
    checkException(fixedFields(vcf) <- as.matrix(values(fixedFields(vcf))), 
        silent=TRUE)

    ## info
    checkTrue(class(values(info(vcf))) == "DataFrame")
    checkException(info(vcf) <- NULL, silent=TRUE)
    v1 <- vcf 
    info(v1) <- DataFrame()
    checkException(validObject(v1), silent=TRUE)

    ## geno
    checkTrue(class(geno(vcf)) == "SimpleList")
    checkException(geno(vcf) <- NULL, silent=TRUE)
    checkIdentical(geno(vcf)$CNQ, geno(vcf)[["CNQ"]])
    checkIdentical(geno(vcf)[["CN"]], assays(vcf)[["CN"]]) 
}

test_VCF_subset <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")

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
    dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                          letters[seq_len(ncol(ss1))])
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
    vcf <- readVcf(fl, "hg19")
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
