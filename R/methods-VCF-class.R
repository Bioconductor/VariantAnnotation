## =========================================================================
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
    if (collapsed)
        class <- "CollapsedVCF"
    else
        class <- "ExpandedVCF"

    new(class, SummarizedExperiment(assays=geno, rowData=rowData,
        colData=colData, exptData=exptData), fixed=fixed, info=info, ...)
}
 
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

### ref 
setMethod("ref", "VCF", 
    function(x) 
{
    slot(x, "fixed")$REF
})

setReplaceMethod("ref", c("VCF", "DNAStringSet"),
    function(x, value)
{
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
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

setReplaceMethod("alt", c("ExpandedVCF", "character"),
    function(x, value)
{
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

### qual 
setMethod("qual", "VCF", 
    function(x) 
{
    slot(x, "fixed")$QUAL
})

setReplaceMethod("qual", c("VCF", "integer"),
    function(x, value)
{
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
    if (!is(value, "GRanges"))
        stop("'value' must be a GRanges")
    idx <- names(mcols(value)) %in% "paramRangeID"
    fixed(x) <- mcols(value)[!idx] 
    slot(x, "rowData") <- value[,idx]
    validObject(x)
    x
})

### We need a 'mcols<-,VCF' distinct from
### 'mcols<-,SummarizedExperiment' because of the
### behavior of 'rowData<-,VCF'. The 'fixed' fields
### are stored in a separate slot and not in
### the meta data columns of the 'rowData' slot.
setReplaceMethod("mcols", "VCF",
    function(x, ..., value)
{
    if (!is(value, "DataFrame"))
        stop("'value' must be a DataFrame")
    idx <- names(value) %in% "paramRangeID"
    fixed(x) <- value[!idx] 
    slot(x, "rowData") <- value[,idx]
    x
})
 
### info 
setMethod("info", "VCF", 
    function(x) 
{
    info <- slot(x, "info")
    if (any(duplicated(rownames(x))))
        return(info)
    if (length(info) != 0L)
        rownames(info) <- rownames(x)
    info
})

setReplaceMethod("info", c("VCF", "DataFrame"),
    function(x, value)
{
    slot(x, "info") <- value
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
    x
})

setReplaceMethod("geno", c("VCF", "character", "matrix"),
    function(x, i, ..., value)
{
    assay(x, i) <- value
    x
})

setReplaceMethod("geno", c("VCF", "numeric", "matrix"),
    function(x, i, ..., value)
{
    assay(x, i) <- value
    x
})

setReplaceMethod("geno", c("VCF", "missing", "matrix"),
    function(x, i, ..., value)
{
    assay(x) <- value
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
    se <- GenomicRanges:::.rbind.SummarizedExperiment(args)
    info <- do.call(rbind, lapply(args, info))
    fixed <- do.call(rbind, lapply(args, fixed))
    new(class(args[[1]]), se, fixed=fixed, info=info)
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
    se <- GenomicRanges:::.cbind.SummarizedExperiment(args)
    info <- GenomicRanges:::.cbind.DataFrame(args, info, "info") 
    if (!.compare(lapply(args, "fixed")))
        stop("data in 'fixed(VCF)' must match.")
    else
        fixed <- fixed(args[[1]]) 
    new(class(args[[1]]), se, fixed=fixed, info=info)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods 
###

setMethod("renameSeqlevels",  c("VCF", "character"),
    function(x, value, ...)
{
    rd <- renameSeqlevels(rowData(x), value)
    rowData(x) <- rd
    x 
})

setMethod("keepSeqlevels",  c("VCF", "character"),
    function(x, value, ...)
{
    rd <- keepSeqlevels(rowData(x), value)
    idx <- as.vector(seqnames(rowData(x))) %in% value
    xsub <- x[idx, ]
    rowData(xsub) <- rd
    xsub 
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
    prettydescr <- function(desc) {
        desc <- as.character(desc)
        wd <- options()[["width"]] - 30
        dwd <- nchar(desc)
        desc <- substr(desc, 1, wd)
        idx <- dwd > wd
        substr(desc[idx], wd - 2, dwd[idx]) <- "..."
        desc
    }
    headerrec <- function(df, lbl, margin="  ") {
        df$Description <- prettydescr(df$Description)
        rownames(df) <- paste(margin, rownames(df))
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
    if (length(hdr <- info(header(object))) > 0) {
        if (nrow(df <- as.data.frame(hdr[colnames(info(object)),])) > 0) {
            cat("info(header(vcf)):\n")
            headerrec(df, "info")
        }
    }
    cat("geno(vcf):\n")
    geno <- geno(object, withDimnames=FALSE)
    printSimpleList(geno, margin=margin) 
    if (nrow(hdr <- geno(header(object))) > 0) {
        if (nrow(df <- as.data.frame(hdr[names(geno),])) > 0) {
            cat("geno(header(vcf)):\n")
            headerrec(df, "geno")
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
