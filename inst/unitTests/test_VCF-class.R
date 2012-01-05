fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")

test_VCF_construction <- function() {
    ## empty-ish
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(VCF()))
    checkTrue(validObject(new("VCF")))
    checkTrue(validObject(VCF(assays=SimpleList(m1))))

    ## substance
    checkTrue(validObject(vcf))
}

test_VCF_getters <- function() {
   checkTrue(class(info(vcf)) == "SimpleList")
   checkTrue(class(info(vcf)[["CIEND"]]) == "array")
   checkIdentical(info(vcf)$BKPTID, info(vcf)[["BKPTID"]])

   checkIdentical(geno(vcf)$CNQ, geno(vcf)[["CNQ"]])
   checkIdentical(geno(vcf)[["CN"]], assays(vcf)[["CN"]]) 
}

#test_VCF_setters <- function()
#{
#}

test_VCF_subset <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")

    ## numeric
    ss1 <- vcf[2:3,]
    checkIdentical(c(2L, ncol(vcf)), dim(ss1))
    ss1 <- vcf[,2]
    checkIdentical(c(nrow(vcf), 1L), dim(ss1))
    checkIdentical(rownames(ss1), names(info(ss1)[[1]]))
    ss1 <- vcf[2:3, 2]
    checkIdentical(c(2L, 1L), dim(ss1))
    checkIdentical(rownames(ss1), names(info(ss1)[[1]]))

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
