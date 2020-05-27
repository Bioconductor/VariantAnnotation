
test_VCF_construction <- function() {
    ## empty
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(VCF()))
    checkTrue(validObject(VCF(collapsed=TRUE)))
    checkTrue(validObject(VCF(collapsed=FALSE)))
    checkTrue(validObject(VCF(geno=SimpleList(m1))))

    ## substance
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    target <- readVcf(fl, genome="hg19")
    current <- VCF(rowRanges=rowRanges(target), colData=colData(target),
                   geno=geno(target), info=info(target),
                   fixed=fixed(target)) 
    checkTrue(validObject(current)) 
    checkIdentical("hg19", unique(genome(rowRanges(target)))) 
}

test_VCF_fixed <- function() {
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")
    vcf1 <- vcf

    ## ref
    checkTrue(class(ref(vcf)) == "DNAStringSet")
    checkException(ref(vcf) <- NULL, silent=TRUE)
    checkException(ref(vcf) <- as.character(ref(vcf)), silent=TRUE)
    ref(vcf1) <- ref(vcf)[1]
    checkIdentical(rep(ref(vcf)[1], nrow(vcf)), ref(vcf1))

    ## alt 
    checkTrue(class(alt(vcf)) == "DNAStringSetList")
    checkException(alt(vcf) <- NULL, silent=TRUE)
    v1 <- vcf
    alt(v1) <- CharacterList(list("A", "A", c("G", "T"), "", c("G", "GTCT"))) 
    checkTrue(validObject(v1))
    alt(vcf1) <- alt(vcf)[1]
    checkIdentical(rep(alt(vcf)[1], nrow(vcf)), alt(vcf1))

    ## qual 
    checkTrue(class(qual(vcf)) == "numeric")
    checkException(qual(vcf) <- NULL, silent=TRUE)
    checkException(qual(vcf) <- as.character(qual(vcf)), silent=TRUE)
    qual(vcf1) <- qual(vcf)[1] # recylcing
    checkIdentical(rep(qual(vcf)[1], nrow(vcf)), qual(vcf1))

    ## filt 
    checkTrue(class(filt(vcf)) == "character")
    checkException(filt(vcf) <- NULL, silent=TRUE)
    checkException(filt(vcf) <- as.list(filt(vcf)), silent=TRUE)
    filt(vcf1) <- filt(vcf)[1] # recylcing
    checkIdentical(rep(filt(vcf)[1], nrow(vcf)), filt(vcf1))

    ## fixed 
    checkTrue(is(fixed(vcf), "DataFrame"))
    checkException(fixed(vcf) <- NULL, silent=TRUE)
    checkException(fixed(vcf) <- as.matrix(fixed(vcf)), 
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
}

test_VCF_rowRanges_info_geno <- function() {
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")
    ## rowRanges
    v1 <- vcf
    rowRanges(v1) <- rowRanges(v1, fixed=FALSE)[5:1]
    checkIdentical(rownames(vcf), rev(rownames(v1)))

    ## info
    checkTrue(is(info(vcf), "DataFrame"))
    checkException(info(vcf) <- NULL, silent=TRUE)
    v1 <- vcf 
    checkException(info(v1) <- DataFrame(), silent=TRUE)

    checkTrue(class(info(vcf)$AF) == "CompressedNumericList")
    AF <- NumericList(0.5, 0.017, c(0.333,0.667), NA, c(NA, NA))
    checkIdentical(info(vcf)$AF, AF) 

    ## geno
    checkTrue(class(geno(vcf)) == "SimpleList")
    checkException(geno(vcf) <- NULL, silent=TRUE)
    checkIdentical(geno(vcf)$CNQ, geno(vcf)[["CNQ"]])
    checkIdentical(geno(vcf)[["CN"]], assays(vcf)[["CN"]])

    checkTrue(is(geno(vcf, "GT"), "matrix"))
    v1 <- vcf
    geno(v1, "GT") <- matrix(NA, nrow=5, ncol=3)
    checkIdentical(unique(as.vector(geno(v1, "GT"))), NA)
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
    checkIdentical(rownames(ss1), rownames(info(ss1)))
    ss1 <- vcf[2:3, 2]
    checkIdentical(c(2L, 1L), dim(ss1))
    checkIdentical(rownames(ss1), rownames(info(ss1)))

    ## character
    ss1 <- vcf 
    dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))], letters[seq_len(ncol(ss1))]) 
    ridx <- c("B", "C") 
    checkIdentical(rowRanges(ss1[ridx,]), rowRanges(ss1)[ridx,]) 
    checkIdentical(rowRanges(ss1["C",]), rowRanges(ss1)["C",,drop=FALSE])
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
    checkIdentical(rowRanges(ss1)[idx,,drop=FALSE], rowRanges(ss2))
    ss2 <- ss1[,idx]
    checkIdentical(colData(ss1)[idx,,drop=FALSE], colData(ss2))

    ## 0 columns
    vcf <- VCF(rowRanges=GRanges("chr1", IRanges(1:10, width=1)))
    checkIdentical(dim(vcf[1:5, ]), c(5L, 0L))
    ## 0 rows 
    vcf <- VCF(colData=DataFrame(samples=1:10))
    checkIdentical(dim(vcf[ ,1:5]), c(0L, 5L))
}

test_VCF_subset_empty_slots <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19", param=ScanVcfParam(info=NA, fixed=NA))
    checkTrue(all(names(mcols(vcf)) %in% c("paramRangeID", "REF"))) 
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
    checkIdentical(rowRanges(vcf)[2:1,], rowRanges(ss1)[1:2,])
    checkIdentical(rowRanges(vcf[-(1:2),]), rowRanges(ss1)[-(1:2),])
    checkIdentical(colData(vcf), colData(ss1))
    checkIdentical(metadata(vcf), metadata(ss1))

    ss1 <- vcf 
    ss1[,1:2] <- ss1[,2:1,drop=FALSE]
    checkIdentical(colData(vcf)[2:1,,drop=FALSE],
                   colData(ss1)[1:2,,drop=FALSE])
    checkIdentical(colData(vcf)[-(1:2),,drop=FALSE],
                   colData(ss1)[-(1:2),,drop=FALSE])
    checkIdentical(rowRanges(vcf), rowRanges(ss1))
    checkIdentical(metadata(vcf), metadata(ss1))
}

test_VCF_seqlevels <- function()
{
    fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, genome="hg19")
    seqlev <- seqlevels(vcf)
    checkIdentical(seqlev, c("1", "2", "3", "4"))

    seqlevels(vcf)[seqlevels(vcf) == "1"] <- "chr1"
    checkIdentical(seqlevels(vcf), c("chr1", "2", "3", "4"))
 
    seqlevels(vcf, pruning.mode="coarse") <- "3"
    checkIdentical(seqlevels(vcf), "3")
}

quiet <- suppressWarnings
test_VCF_cbind <- function()
## requires matching ranges
{
    ## empty
    vcf <- VCF()
    empty <- cbind(vcf, vcf)
    checkEquals(dim(empty), c(0, 0))

    ## different ranges 
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf1 <- readVcf(fl, genome="hg19")
    vcf2 <- vcf1[2:4]
    rownames(vcf2) <- month.name[seq_len(nrow(vcf2))]
    checkException(quiet(cbind(vcf1, vcf2)), silent=TRUE)

    ## same ranges
    vcf2 <- vcf1[,1]
    #colnames(vcf2) <- month.name[seq_len(ncol(vcf2))]
    res <- cbind(vcf1, vcf2)
    checkTrue(nrow(res) == 5)
    checkTrue(ncol(res) == 4)
    ## info 
    checkTrue(ncol(info(res)) == 6)
    info2 <- info(vcf2)
    info2$H2 <- !info2$H2
    info(vcf2) <- info2
    checkException(cbind(vcf1, vcf2), silent=TRUE)
    info(vcf2) <- info(vcf2)[,1:2] 
    res <- cbind(vcf1, vcf2)
    checkTrue(ncol(info(res)) == 6) 
    info(vcf2) <- DataFrame("noMatch"=1:5)
    res <- cbind(vcf1, vcf2)
    checkTrue(ncol(info(res)) == 7) 
    ## fixed 
    vcf2 <- vcf1[,1]
    res <- cbind(vcf1, vcf2)
    checkTrue(ncol(fixed(res)) == 4)
    fixed2 <- fixed(vcf2)
    fixed2$QUAL <- c(rep(1, 5)) 
    fixed(vcf2) <- fixed2
    checkException(cbind(vcf1, vcf2), silent=TRUE)
}

test_VCF_rbind <- function()
## requires matching samples 
{
    ## empty
    vcf0 <- VCF()
    empty <- rbind(vcf0, vcf0)
    checkEquals(dim(empty), c(0, 0))

    ## different sample 
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf1 <- readVcf(fl, genome="hg19")
    vcf2 <- vcf1[,2]
    colnames(vcf2) <- month.name[seq_len(ncol(vcf2))]
    checkException(quiet(rbind(vcf1, vcf2)), silent=TRUE)

    ## same samples 
    vcf2 <- vcf1[1:3,]
    res <- rbind(vcf1, vcf2)
    checkTrue(nrow(res) == 8)
    checkTrue(ncol(res) == 3)
    ## info 
    checkTrue(ncol(info(res)) == 6)
    info2 <- info(vcf2)
    info2$H2 <- !info2$H2
    info(vcf2) <- info2
    res <- rbind(vcf1, vcf2)
    checkTrue(ncol(info(res)) == 6)
    ## fixed 
    vcf2 <- vcf1[1:3,]
    res <- rbind(vcf1, vcf2)
    checkTrue(ncol(fixed(res)) == 4)
    fixed2 <- fixed(vcf2)
    fixed2$QUAL <- c(rep(1, nrow(vcf2))) 
    fixed(vcf2) <- fixed2
    res <- rbind(vcf1, vcf2)
    checkTrue(ncol(fixed(res)) == 4)
}
