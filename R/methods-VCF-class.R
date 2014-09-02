### =========================================================================
### VCF class methods 
### =========================================================================


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###

VCF <-
    function(rowData=GRanges(), colData=DataFrame(), 
             exptData=SimpleList(header=VCFHeader()), 
             fixed=DataFrame(),
             info=DataFrame(row.names=seq_along(rowData)), 
             geno=SimpleList(),
             ..., collapsed=TRUE, verbose=FALSE)
{
    rownames(info) <- rownames(fixed) <- NULL
    if (collapsed) {
        class <- "CollapsedVCF"
        if (!length(fixed)) {
            fixed=DataFrame(
                REF=DNAStringSet(rep(".", length(rowData))), 
                ALT=DNAStringSetList(as.list(rep(".", length(rowData)))),
                QUAL=rep(NA_real_, length(rowData)), 
                FILTER=rep(NA_character_, length(rowData)))
        }
    } else {
        class <- "ExpandedVCF"
        if (!length(fixed)) {
            fixed=DataFrame(
                REF=DNAStringSet(rep("", length(rowData))), 
                ALT=DNAStringSet(rep("", length(rowData))),
                QUAL=rep(NA_real_, length(rowData)), 
                FILTER=rep(NA_character_, length(rowData)))
        }
    }

    new(class, SummarizedExperiment(assays=geno, rowData=rowData,
        colData=colData, exptData=exptData), fixed=fixed, info=info, ...)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

.checkLength <- function (x, len) 
{
    if (len != length(slot(x, "rowData"))) 
        stop("length(value) must equal length(rowData(x))")
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

### rowData
setMethod("rowData", "VCF", 
    function(x) 
{
    gr <- slot(x, "rowData") 
    if (length(slot(x, "fixed")) != 0L)
        mcols(gr) <- append(mcols(gr), slot(x, "fixed"))
    gr
})

setReplaceMethod("rowData", c("VCF", "GRanges"),
    function(x, value)
{
    fixed_idx <- names(mcols(value)) %in% c("REF", "ALT", "QUAL", "FILTER")
    slot(x, "fixed") <- mcols(value)[fixed_idx] 
    slot(x, "rowData") <- value[,!fixed_idx]
    validObject(x)
    x
})

### Must define VCF methods for 'mcols<-' and 'dimnames<-'
### instead of inheriting from SummarizedExperiment.
### The 'fixed' fields are stored in a separtate slot
### and not as metadata columns of 'rowData'.

setReplaceMethod("mcols", c("VCF", "DataFrame"),
    function(x, ..., value)
{
    idx <- names(value) %in% "paramRangeID"
    fixed(x) <- value[!idx] 
    slot(x, "rowData") <- value[,idx]
    x
})

setReplaceMethod("dimnames", c("VCF", "list"),
    function(x, value)
{
    rowData <- slot(x,"rowData")
    names(rowData) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    GenomicRanges:::clone(x, rowData=rowData, colData=colData)
})
 
### info 
setMethod("info", "VCF", 
    function(x, ..., row.names = TRUE) 
{
    info <- slot(x, "info")
    if (row.names) {
        if (length(info) != 0L)
            if (!sum(duplicated(rownames(x))))
                rownames(info) <- rownames(x)
    }
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
    strand(rowData(x))
})

setReplaceMethod("strand", "VCF",
    function(x, ..., value)
{
    strand(rowData(x)) <- value
    x
})

### header
setMethod("header", "VCF",
    function(x)
{ 
    exptData(x)$header
})

setReplaceMethod("header", c("VCF", "VCFHeader"),
    function(x, value)
{
    slot(x, "exptData")$header <- value
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
        i <- GenomicRanges:::.SummarizedExperiment.charbound(i, rownames(x), msg)
    }

    ii <- ff <- NULL
    if (!missing(i)) {
        if (length(slot(x, "info")) != 0L)
            ii <- i
        if (length(slot(x, "fixed")) != 0L)
            ff <- i
    }
    if (missing(i) && missing(j)) {
        x
    } else if (missing(i)) {
        callNextMethod(x, , j, ...)
    } else if (missing(j)) {
        callNextMethod(x, i, , rowData=slot(x, "rowData")[i,,drop=FALSE],
                       info=slot(x, "info")[ii,,drop=FALSE],
                       fixed=slot(x, "fixed")[ff,,drop=FALSE], ...)
    } else {
        callNextMethod(x, i, j, rowData=slot(x, "rowData")[i,,drop=FALSE],
                       info=slot(x, "info")[ii,,drop=FALSE],
                       fixed=slot(x, "fixed")[ff,,drop=FALSE], ...)
    }
})

setReplaceMethod("[",
    c("VCF", "ANY", "ANY", "VCF"),
    function(x, i, j, ..., value)
{
    if (!missing(i) && is.character(i)) {
        msg <- "<VCF>[i,]<- index out of bounds: %s"
        i <- GenomicRanges:::.SummarizedExperiment.charbound(i, rownames(x), msg)
    }

    if (missing(i) && missing(j)) {
        x
    } else if (missing(i)) {
        callNextMethod(x, , j, ..., value=value)
    } else if (missing(j)) {
        callNextMethod(x, i, , rowData=local({
            rd <- slot(x, "rowData")
            rd[i,] <- slot(value, "rowData")
            names(rd)[i] <- names(slot(value, "rowData"))
            rd
        }), info=local({
            ii <- slot(x, "info")
            ii[i,] <- slot(value, "info")
            ii 
        }), fixed=local({
            ff <- slot(x, "fixed") 
            ff[i,] <- slot(value, "fixed") 
            ff
        }), ..., value=value)
    } else {
        callNextMethod(x, i, j, rowData=local({
            rd <- slot(x, "rowData")
            rd[i,] <- slot(value, "rowData")
            names(rd)[i] <- names(slot(value, "rowData"))
            rd
        }), info=local({
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

.compare <- GenomicRanges:::.compare
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

    rowData <- do.call(c, lapply(args,
        function(i) slot(i, "rowData")))
    colData <- GenomicRanges:::.cbind.DataFrame(args, colData, "colData")
    assays <- GenomicRanges:::.bind.arrays(args, rbind, "assays")
    exptData <- do.call(c, lapply(args, exptData))
    info <- do.call(rbind, lapply(args, info))
    fixed <- do.call(rbind, lapply(args, fixed))
    new(class(args[[1]]), assays=assays, rowData=rowData,
        colData=colData, exptData=exptData, fixed=fixed, info=info)
})

.renameSamples <- function(args)
{
    Map(function(x, i) {
        idx <- names(colData(x)) == "Samples"
        names(colData(x))[idx] <- paste0("Samples.", i)
        x
    }, x=args, i=seq_along(args)) 
}

## Appropriate for objects with same ranges and different samples.
setMethod("cbind", "VCF",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    if (!.compare(lapply(args, class)))
        stop("'...' objects must be of the same VCF class")


    if (!.compare(lapply(args, rowData), TRUE))
        stop("'...' object ranges (rows) are not compatible")
    rowData <- rowData(args[[1]])
    mcols(rowData) <- GenomicRanges:::.cbind.DataFrame(args, mcols, "mcols")
    info <- GenomicRanges:::.cbind.DataFrame(args, info, "info") 
    colData <- do.call(rbind, lapply(args, colData))
    assays <- GenomicRanges:::.bind.arrays(args, cbind, "assays")
    exptData <- do.call(c, lapply(args, exptData))
    if (!.compare(lapply(args, "fixed")))
        stop("data in 'fixed(VCF)' must match.")
    else
        fixed <- fixed(args[[1]]) 
    new(class(args[[1]]), assays=assays, rowData=rowData,
        colData=colData, exptData=exptData, fixed=fixed, info=info)


})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show methods
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
    cat("rowData(vcf):\n")
    printSmallGRanges(rowData(object), margin=margin)
    cat("info(vcf):\n")
    printSmallDataTable(info(object), margin=margin) 
    if (length(header(object))) {
        if (length(hdr <- info(header(object)))) {
            nms <- intersect(colnames(info(object)), rownames(hdr))
            diff <- setdiff(colnames(info(object)), rownames(hdr))
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
### updateObject 
###

### Creates a CollapsedVCF from a VCF.
setMethod("updateObject", "VCF",
    function(object, ..., verbose=FALSE)
    {
        if (verbose) 
            message("updateObject(object = 'VCF')") 
        VCF(rowData=rowData(object), colData=colData(object), exptData=exptData(object), 
            info=mcols(info(object))[-1], fixed=mcols(fixed(object))[-1], geno=geno(object))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### restrictToSNV 
###

### Creates a subset VCF with SNVs only. 'x' can be a CollapsedVCF 
### or ExpandedVCF. Data in 'alt(x)' must be DNAStringSet or DNAStringSetList.
restrictToSNV <- function(x, ...)
{
    .Deprecated("isSNV")
}
