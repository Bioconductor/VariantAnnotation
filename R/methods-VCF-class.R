### =========================================================================
### VCF class methods 
### =========================================================================


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###

VCF <-
    function(rowRanges=GRanges(), colData=DataFrame(), 
             exptData=list(header=VCFHeader()), 
             fixed=DataFrame(),
             info=DataFrame(row.names=seq_along(rowRanges)), 
             geno=SimpleList(),
             ..., collapsed=TRUE, verbose=FALSE)
{
    rownames(info) <- rownames(fixed) <- NULL
    if (collapsed) {
        class <- "CollapsedVCF"
        if (!length(fixed)) {
            fixed=DataFrame(
                REF=DNAStringSet(rep(".", length(rowRanges))), 
                ALT=DNAStringSetList(as.list(rep(".", length(rowRanges)))),
                QUAL=rep(NA_real_, length(rowRanges)), 
                FILTER=rep(NA_character_, length(rowRanges)))
        }
    } else {
        class <- "ExpandedVCF"
        if (!length(fixed)) {
            fixed=DataFrame(
                REF=DNAStringSet(rep("", length(rowRanges))), 
                ALT=DNAStringSet(rep("", length(rowRanges))),
                QUAL=rep(NA_real_, length(rowRanges)),
                FILTER=rep(NA_character_, length(rowRanges)))
        }
    }

    new(class, SummarizedExperiment(assays=geno, rowRanges=rowRanges,
        colData=colData, metadata=as.list(exptData)), fixed=fixed, info=info,
        ...)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

.checkLength <- function (x, len) 
{
    if (len != length(slot(x, "rowRanges"))) 
        stop("length(value) must equal length(rowRanges(x))")
}

### ref 
setMethod("ref", "VCF", 
    function(x) 
{
    slot(x, "fixed")$REF
})

setReplaceMethod("ref", c("VCF", "DNAStringSet"),
    function(x, value)
{
    .checkLength(x, length(value))
    slot(x, "fixed")$REF <- value
    x
})

### alt 
setMethod("alt", "VCF", 
    function(x) 
{
    slot(x, "fixed")$ALT
})

setReplaceMethod("alt", c("CollapsedVCF", "CharacterList"),
    function(x, value)
{
    .checkLength(x, length(value))
    slot(x, "fixed")$ALT <- value
    x
})

setReplaceMethod("alt", c("ExpandedVCF", "character"),
    function(x, value)
{
    .checkLength(x, length(value))
    slot(x, "fixed")$ALT <- value
    x
})

### qual 
setMethod("qual", "VCF", 
    function(x) 
{
    slot(x, "fixed")$QUAL
})

setReplaceMethod("qual", c("VCF", "numeric"),
    function(x, value)
{
    .checkLength(x, length(value))
    slot(x, "fixed")$QUAL <- value
    x
})

### filt
setMethod("filt", "VCF", 
    function(x) 
{
    slot(x, "fixed")$FILTER
})

setReplaceMethod("filt", c("VCF", "character"),
    function(x, value)
{
    .checkLength(x, length(value))
    slot(x, "fixed")$FILTER <- value
    x
})

### fixed 
setMethod("fixed", "VCF", 
    function(x) 
{
    slot(x, "fixed")
})

setReplaceMethod("fixed", c("VCF", "DataFrame"),
    function(x, value)
{
    slot(x, "fixed") <- value
    validObject(x)
    x
})

### rowRanges
setMethod("rowRanges", "VCF", 
    function(x, ..., fixed = TRUE) 
{
    gr <- slot(x, "rowRanges")
    if (fixed)
        if (length(slot(x, "fixed")) != 0L)
            mcols(gr) <- append(mcols(gr), slot(x, "fixed"))
    gr
})

setReplaceMethod("rowRanges", c("VCF", "GRanges"),
    function(x, value)
{
    fixed_idx <- names(mcols(value)) %in% c("REF", "ALT", "QUAL", "FILTER")
    slot(x, "fixed") <- mcols(value)[fixed_idx] 
    slot(x, "rowRanges") <- value[,!fixed_idx]
    validObject(x)
    x
})

### Must define VCF methods for 'mcols<-' and 'dimnames<-' instead of 
### inheriting from RangedSummarizedExperiment. 'fixed' fields are stored in 
### a separate slot and not as metadata columns of 'rowRanges'.

setReplaceMethod("mcols", c("VCF", "DataFrame"),
    function(x, ..., value)
{
    idx <- names(value) %in% "paramRangeID"
    fixed(x) <- value[!idx] 
    slot(x, "rowRanges") <- value[,idx]
    x
})

setReplaceMethod("dimnames", c("VCF", "list"),
    function(x, value)
{
    rowRanges <- slot(x, "rowRanges")
    names(rowRanges) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    BiocGenerics:::replaceSlots(x, rowRanges=rowRanges, colData=colData)
})
 
### info 
setMethod("info", "VCF", 
    function(x, ..., row.names = TRUE) 
{
    info <- slot(x, "info")
    if (row.names && !anyDuplicated(rownames(x)))
        rownames(info) <- rownames(x)
    info
})

setReplaceMethod("info", c("VCF", "DataFrame"),
    function(x, value)
{
    slot(x, "info") <- value
    .valid.VCFHeadervsVCF.fields(x, info)
    validObject(x)
    x
})

### geno
### 'assays' extracts full list
### 'assay' extracts individual list elements

setMethod("geno", "VCF",
    function(x, ..., withDimnames = TRUE)
{
    assays(x, ..., withDimnames=withDimnames) 
})

setMethod("geno", c("VCF", "numeric"),
    function(x, i, ..., withDimnames = TRUE)
{
    assay(x, i, ...) 
})

setMethod("geno", c("VCF", "character"),
    function(x, i, ..., withDimnames = TRUE)
{
    assay(x, i, ...) 
})

setReplaceMethod("geno", c("VCF", "missing", "SimpleList"),
    function(x, i, ..., value)
{
    assays(x) <- value
    .valid.VCFHeadervsVCF.fields(x, geno)
    x
})

setReplaceMethod("geno", c("VCF", "character", "matrix"),
    function(x, i, ..., value)
{
    assay(x, i) <- value
    .valid.VCFHeadervsVCF.fields(x, geno)
    x
})

setReplaceMethod("geno", c("VCF", "numeric", "matrix"),
    function(x, i, ..., value)
{
    assay(x, i) <- value
    .valid.VCFHeadervsVCF.fields(x, geno)
    x
})

setReplaceMethod("geno", c("VCF", "missing", "matrix"),
    function(x, i, ..., value)
{
    assay(x) <- value
    .valid.VCFHeadervsVCF.fields(x, geno)
    x
})


### strand
setMethod("strand", "VCF",
    function(x, ...)
{
    strand(rowRanges(x))
})

setReplaceMethod("strand", "VCF",
    function(x, ..., value)
{
    strand(rowRanges(x)) <- value
    x
})

### header
setMethod("header", "VCF",
    function(x)
{ 
    metadata(x)$header
})

setReplaceMethod("header", c("VCF", "VCFHeader"),
    function(x, value)
{
    slot(x, "metadata")$header <- value
    .valid.VCFHeadervsVCF(x)
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining
###

setMethod("[", c("VCF", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,VCF,ANY,ANY-method'")

    if (!missing(i) && is.character(i)) {
        msg <- "<VCF>[i,] index out of bounds: %s"
        i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i, rownames(x), msg)
    }

    if (missing(i) && missing(j)) {
        x
    } else if (missing(i)) {
        callNextMethod(x, , j, ...)
    } else if (missing(j)) {
        callNextMethod(x, i, ,
                       info=slot(x, "info")[i,,drop=FALSE],
                       fixed=slot(x, "fixed")[i,,drop=FALSE], ...)
    } else {
        callNextMethod(x, i, j,
                       info=slot(x, "info")[i,,drop=FALSE],
                       fixed=slot(x, "fixed")[i,,drop=FALSE], ...)
    }
})

setReplaceMethod("[",
    c("VCF", "ANY", "ANY", "VCF"),
    function(x, i, j, ..., value)
{
    if (!missing(i) && is.character(i)) {
        msg <- "<VCF>[i,]<- index out of bounds: %s"
        i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i, rownames(x), msg)
    }

    if (missing(i) && missing(j)) {
        x
    } else if (missing(i)) {
        callNextMethod(x, , j, ..., value=value)
    } else if (missing(j)) {
        callNextMethod(x, i, ,
            info=local({
            ii <- slot(x, "info")
            ii[i,] <- slot(value, "info")
            ii 
        }), fixed=local({
            ff <- slot(x, "fixed") 
            ff[i,] <- slot(value, "fixed") 
            ff
        }), ..., value=value)
    } else {
        callNextMethod(x, i, j,
            info=local({
            ii <- slot(x, "info")
            ii[i,] <- slot(value, "info")
            ii 
        }), fixed=local({
            ff <- slot(x, "fixed") 
            ff[i,] <- slot(value, "fixed") 
            ff
        }), ..., value=value)
    }
})

.compare <- SummarizedExperiment:::.compare
## Appropriate for objects with different ranges and same samples.
setMethod("rbind", "VCF",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    if (!.compare(lapply(args, class)))
        stop("'...' objects must be of the same VCF class")
    args <- .renameSamples(args)
    if (!.compare(lapply(args, colnames)))
            stop("'...' objects must have the same colnames")
    if (!.compare(lapply(args, ncol)))
            stop("'...' objects must have the same number of samples")

    fixed <- do.call(rbind, lapply(args, fixed))
    info <- do.call(rbind, lapply(args, info))
    ## rowRanges should not contain REF, ALT, QUAL or FILTER
    rowRanges <- do.call(c, lapply(args, slot, "rowRanges"))
    colData <- SummarizedExperiment:::.cbind.DataFrame(args, colData, "colData")
    assays <- do.call(rbind, lapply(args, slot, "assays"))
    elementMetadata <- do.call(rbind, lapply(args, slot, "elementMetadata"))
    metadata <- do.call(c, lapply(args, metadata))

    BiocGenerics:::replaceSlots(args[[1L]],
        fixed=fixed, info=info,
        rowRanges=rowRanges, colData=colData, assays=assays,
        elementMetadata=elementMetadata, metadata=metadata)
})

.renameSamples <- function(args)
{
    Map(function(x, i) {
        idx <- names(colData(x)) == "Samples"
        names(colData(x))[idx] <- paste0("Samples.", i)
        x
    }, x=args, i=seq_along(args)) 
}

## FIXME: combine header info
## Appropriate for objects with same ranges and different samples.
setMethod("cbind", "VCF",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    if (!.compare(lapply(args, class)))
        stop("'...' objects must be of the same VCF class")
    if (!.compare(lapply(args, rowRanges), TRUE))
        stop("'...' object ranges (rows) are not compatible")
    if (!.compare(lapply(args, "fixed")))
        stop("data in 'fixed(VCF)' must match.")

    fixed <- fixed(args[[1L]])
    info <- SummarizedExperiment:::.cbind.DataFrame(args, info, "info") 
    rowRanges <- rowRanges(args[[1L]])
    mcols(rowRanges) <- SummarizedExperiment:::.cbind.DataFrame(args, mcols, "mcols")
    colData <- do.call(rbind, lapply(args, colData))
    assays <- do.call(cbind, lapply(args, slot, "assays"))
    metadata <- do.call(c, lapply(args, metadata))

    BiocGenerics:::replaceSlots(args[[1L]],
        fixed=fixed, info=info,
        rowRanges=rowRanges, colData=colData, assays=assays,
        metadata=metadata)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion 
###

SnpMatrixToVCF <- function(from, seqSource)
{
    if (!all(names(from) %in% c("genotypes", "fam", "map")))
        stop("'from' must be a list with named elements 'genotypes', ",
             "'fam' and 'map'")

    ## match seqlevels style to seqSource
    map <- from$map
    chr <- as.character(map$chromosome)
    seqlevelsStyle(chr) <- seqlevelsStyle(seqSource)
    uniqueChr <- unique(chr)
    if (any(invalid <- !uniqueChr %in% seqlevels(seqSource)))
        stop("seqlevels not found in 'seqSource': ", 
             paste(uniqueChr[invalid], collapse=", "))

    ## ref
    gr <- GRanges(chr, IRanges(map$position, width=1))
    ref <- getSeq(seqSource, gr)

    ## alt 
    refIsAlt1 <- ref == map$allele.1
    refIsAlt2 <- ref == map$allele.2
    alt <- DNAStringSetList(as.list(map$allele.1))
    alt[refIsAlt1] <- DNAStringSetList(as.list(map$allele.2[refIsAlt1]))
    if (any(refNotFound <- refIsAlt1 & refIsAlt2)) {
        df <- data.frame(rbind(map$allele.1[refNotFound], 
                               map$allele.2[refNotFound]))
        alt[refNotFound] <- DNAStringSetList(as.list(df))
    }

    ## genotypes 
    GT <- as.raw(from$genotypes)
    nrowGT <- nrow(from$genotypes)
    newGT <- c("./.", "0/0", "0/1", "1/1") 
    GT <- newGT[as.integer(GT) + 1L]
    GT <- matrix(GT, byrow=TRUE, ncol=nrowGT)

    vcf <- VCF(rowRanges=gr, fixed=DataFrame(REF=ref, ALT=alt), 
               colData=DataFrame(Samples=seq_len(nrowGT)), 
               geno=SimpleList(GT=GT))
    dimnames(vcf) <- list(colnames(from$genotypes), rownames(from$genotypes))
    vcf
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show methods
###

setMethod(show, "VCF",
    function(object)
{
    paste0("This object is no longer valid. Please use updateObject() to ",
           "create a CollapsedVCF instance.")
})

### methods for CollapsedVCF and ExapandedVCF
.showVCFSubclass <- function(object)
{
    prettydescr <- function(desc, wd0=30L) {
        desc <- as.character(desc)
        wd <- options()[["width"]] - wd0
        dwd <- nchar(desc)
        desc <- substr(desc, 1, wd)
        idx <- dwd > wd
        substr(desc[idx], wd - 2, dwd[idx]) <- "..."
        desc
    }
    headerrec <- function(df, margin="   ") {
        wd <- max(nchar(rownames(df))) + nchar("Number") +
            max(nchar(df$Type)) + 7L
        df$Description <- prettydescr(df$Description, wd)
        rownames(df) <- paste0(margin, rownames(df))
        print(df, right=FALSE)
    }
    printSmallGRanges <- function(x, margin)
    {
        lx <- length(x)
        nc <- ncol(mcols(x))
        nms <- names(mcols(x))
        cat(margin, class(x), " with ",
            nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
            ": ", paste0(nms, collapse=", "), "\n", sep="")
    }
    printSmallDataTable <- function(x, margin)
    {
        nc <- ncol(x)
        nms <- names(x)
        cat(margin, class(x), " with ",
            nc, ifelse(nc == 1, " column", " columns"),
            ": ", prettydescr(paste0(nms, collapse=", ")), "\n", sep="")
    }
    printSimpleList <- function(x, margin)
    {
        lo <- length(x)
        nms <- names(x)
        cat(margin, class(x), " of length ", lo, 
            ": ", prettydescr(paste0(nms, collapse=", ")), "\n", sep="")
    }
    margin <- "  "
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat("rowRanges(vcf):\n")
    printSmallGRanges(rowRanges(object), margin=margin)
    cat("info(vcf):\n")
    printSmallDataTable(info(object, row.names=FALSE), margin=margin) 
    if (length(header(object))) {
        if (length(hdr <- info(header(object)))) {
            nms <- intersect(colnames(info(object, row.names=FALSE)),
                             rownames(hdr))
            diff <- setdiff(colnames(info(object, row.names=FALSE)),
                            rownames(hdr))
            if (length(nms)) {
                cat("info(header(vcf)):\n")
                headerrec(as.data.frame(hdr[nms,]))
            }
            if (length(diff))
                cat("  Fields with no header:", 
                    paste(diff, collapse=","), "\n")
        }
    }
    cat("geno(vcf):\n")
    geno <- geno(object, withDimnames=FALSE)
    printSimpleList(geno, margin=margin) 
    if (length(header(object))) {
        if (length(hdr <- geno(header(object)))) {
            nms <- intersect(names(geno(object)), rownames(hdr))
            diff <- setdiff(names(geno(object)), rownames(hdr))
            if (length(nms)) {
                cat("geno(header(vcf)):\n")
                headerrec(as.data.frame(hdr[nms,]))
            }
            if (length(diff))
                cat("  Fields with no header:", 
                    paste(diff, collapse=","), "\n")
        }
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### restrictToSNV() 
###

### Creates a subset VCF with SNVs only. 'x' can be a CollapsedVCF 
### or ExpandedVCF. Data in 'alt(x)' must be DNAStringSet or DNAStringSetList.
restrictToSNV <- function(x, ...)
{
    .Defunct("isSNV")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### genotypeCodesToNucleotides()
###

genotypeCodesToNucleotides <- function(vcf, ...) 
{
    GT <- geno(vcf, withDimnames=FALSE)$GT
    if (is.null(GT <- geno(vcf, withDimnames=FALSE)$GT))
        stop("no 'GT' data found in geno(vcf)")

    if (is(alt(vcf), "CharacterList"))
        stop("'ALT' must be a DNAStringSetList")
    ALT <- as.list(splitAsList(as.character(alt(vcf)@unlistData),
                   togroup(PartitioningByWidth(alt(vcf))))) 
    REF <- as.character(ref(vcf))
    geno(vcf)$GT <- .geno2geno(NULL, ALT, REF, GT)
    vcf
}

