## =========================================================================
### VCF class methods 
### =========================================================================


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###

VCF <-
    function(rowData=GRanges(), colData=DataFrame(), exptData=SimpleList(), 
             fixed=DataFrame(), info=DataFrame(), geno=SimpleList(),
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
setMethod("geno", "VCF",
    function(x, ..., withDimnames = TRUE)
{
    assays(x, ..., withDimnames=withDimnames) 
})

setReplaceMethod("geno", c("VCF", "character", "matrix"),
    function(x, i, ..., value)
{
    assays(x)[[i]] <- value
    x
})

setReplaceMethod("geno", c("VCF", "numeric", "matrix"),
    function(x, i, ..., value)
{
    assays(x)[[i]] <- value
    x
})

setReplaceMethod("geno", c("VCF", "missing", "SimpleList"),
    function(x, i, ..., value)
{
    assays(x) <- value
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
.cbind.DataFrame <- GenomicRanges:::.cbind.DataFrame
.rbind.DataFrame <- GenomicRanges:::.rbind.DataFrame
.cbind.SummarizedExperiment <- GenomicRanges:::.cbind.SummarizedExperiment
.rbind.SummarizedExperiment <- GenomicRanges:::.rbind.SummarizedExperiment
## Appropriate for objects with different ranges and same samples.
setMethod("rbind", "VCF",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    if (!.compare(lapply(args, class)))
        stop("'...' objects must be of the same VCF class")
    se <- .rbind.SummarizedExperiment(args)
    info <- .rbind.DataFrame(args, info, "info") 
    fixed <- .rbind.DataFrame(args, fixed, "fixed") 
    new(class(args[[1]]), se, fixed=fixed, info=info)
})

## Appropriate for objects with same ranges and different samples.
setMethod("cbind", "VCF",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    if (!.compare(lapply(args, class)))
        stop("'...' objects must be of the same VCF class")
    se <- .cbind.SummarizedExperiment(args)
    info <- .cbind.DataFrame(args, info, "info") 
    fixed <- .cbind.DataFrame(args, fixed, "fixed") 
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
### show
###

setMethod(show, "VCF",
    function(object)
{
    paste0("This object is no longer valid. Please use updateObject() to ",
           "create a CollapsedVCF instance.")
})

### show method for CollapsedVCF and ExapandedVCF
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
            ": ", paste0(nms, collapse=", "), "\n", sep = "")
    }
    printSimpleList <- function(x, margin)
    {
        lo <- length(x)
        nms <- names(x)
        cat(margin, class(x), " of length ", lo, 
            ": ", paste0(nms, collapse=", "), "\n", sep = "")
    }
    margin <- "  "
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat("rowData(vcf):\n")
    printSmallGRanges(rowData(object), margin=margin)
    cat("info(vcf):\n")
    printSmallDataTable(info(object), margin=margin) 
    if (!is.null(header(object))) {
    cat("info(header(vcf)):\n")
        info <- info(header(object))[colnames(info(object)),]
        if (nrow(info) > 0)
            headerrec(as.data.frame(info), "info")
    }
    cat("geno(vcf):\n")
    printSimpleList(geno(object), margin=margin) 
    if (!is.null(header(object))) {
    cat("geno(header(vcf)):\n")
        geno <- as.data.frame(geno(exptData(object)$header))
        headerrec(geno, "geno")
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
