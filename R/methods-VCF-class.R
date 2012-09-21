### =========================================================================
### VCF class methods 
### =========================================================================


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Constructor 
##

VCF <-
    function(rowData=GRanges(), colData=DataFrame(), exptData=SimpleList(), 
             fixed=DataFrame(row.names=1:length(rowData)), 
             info=DataFrame(row.names=1:length(rowData)), geno=SimpleList(),
             ..., verbose=FALSE)
{
    rownames(info) <- rownames(fixed) <- NULL
    new("VCF", SummarizedExperiment(assays=geno, rowData=rowData,
        colData=colData, exptData=exptData), fixed=fixed, info=info, ...)
}
 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Getters and Setters
##

## ref 
setMethod("ref", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    if (!is.null(slot(x, "fixed")$REF))
        values(gr) <- append(values(gr), 
            DataFrame(REF=slot(x, "fixed")$REF))
    gr
})

setReplaceMethod("ref", c("VCF", "DNAStringSet"),
    function(x, value)
{
    slot(x, "fixed")$REF <- value
    x
})

## alt 
setMethod("alt", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    if (!is.null(slot(x, "fixed")$ALT))
        values(gr) <- append(values(gr),
            DataFrame(ALT=slot(x, "fixed")$ALT))
    gr
})

setReplaceMethod("alt", c("VCF", "CharacterList"),
    function(x, value)
{
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

setReplaceMethod("alt", c("VCF", "DNAStringSetList"),
    function(x, value)
{
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

## qual 
setMethod("qual", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    if (length(slot(x, "fixed")$QUAL) != 0L)
        values(gr) <- append(values(gr), 
            DataFrame(QUAL=slot(x, "fixed")$QUAL))
    gr
})

setReplaceMethod("qual", c("VCF", "integer"),
    function(x, value)
{
    slot(x, "fixed")$QUAL <- value
    x
})

## filt
setMethod("filt", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    if (length(slot(x, "fixed")$FILTER) != 0L)
        values(gr) <- append(values(gr), 
            DataFrame(FILTER=slot(x, "fixed")$FILTER))
    gr
})

setReplaceMethod("filt", c("VCF", "character"),
    function(x, value)
{
    slot(x, "fixed")$FILTER <- value
    x
})

## fixed 
setMethod("fixed", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    if (length(slot(x, "fixed")) != 0L)
        values(gr) <- append(values(gr), slot(x, "fixed"))
    gr
})

setReplaceMethod("fixed", c("VCF", "DataFrame"),
    function(x, value)
{
    slot(x, "fixed") <- value
    validObject(x)
    x
})

## info 
setMethod("info", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    if (length(slot(x, "info")) != 0L)
        values(gr) <- append(values(gr), slot(x, "info"))
    gr
})

setReplaceMethod("info", c("VCF", "DataFrame"),
    function(x, value)
{
    slot(x, "info") <- value
    validObject(x)
    x
})

## geno
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

## strand
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


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Subsetting 
##

setMethod("[", c("VCF", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,VCF,ANY,ANY-method'")

    if (!missing(i) && is.character(i)) {
        msg <- "<VCF>[i,] index out of bounds: %s"
        i <- GenomicRanges:::.SummarizedExperiment.charbound(i, rownames(x), msg)
    }

    if (missing(i) && missing(j)) {
        x
    } else if (missing(i)) {
        callNextMethod(x, , j, ...)
    } else if (missing(j)) {
        callNextMethod(x, i, , info=slot(x, "info")[i,,drop=FALSE],
                       fixed=slot(x, "fixed")[i,,drop=FALSE], ...)
    } else {
        callNextMethod(x, i, j, info=slot(x, "info")[i,,drop=FALSE],
                       fixed=slot(x, "fixed")[i,,drop=FALSE], ...)
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
        callNextMethod(x, i, , info=local({
            ii <- slot(x, "info")
            ii[i,] <- slot(value, "info")
            ii 
        }), fixed=local({
            ff <- slot(x, "fixed") 
            ff[i,] <- slot(value, "fixed") 
            ff
        }), ..., value=value)
    } else {
        callNextMethod(x, i, j, info=local({
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

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Other methods 
##

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


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Show
##

setMethod(show, "VCF",
    function(object)
{
    selectSome <- IRanges:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat("genome:", unique(genome(rowData(object))), "\n")

    expt <- names(exptData(object))
    if (is.null(expt))
        expt <- character(length(exptData(object)))
    scat("exptData(%d): %s\n", expt)

    fixed <- names(slot(object, "fixed")) 
    if (is.null(fixed))
        fixed <- character(ncol(values(fixed(object))))
    scat("fixed(%d): %s\n", fixed)

    info <- names(slot(object, "info")) 
    if (is.null(info))
        info <- character(ncol(values(info(object))))
    scat("info(%d): %s\n", info)

    nms <- names(geno(object, withDimnames=FALSE))
    if (is.null(nms))
        nms <- character(length(geno(object, withDimnames=FALSE)))
    scat("geno(%d): %s\n", nms)

    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]]) scat("rownames(%d): %s\n", dimnames[[1]])
    else scat("rownames: NULL\n")
    scat("rowData values names(%d): %s\n",
         names(values(rowData(object))))
    if (dlen[[2]]) scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")
    scat("colData names(%d): %s\n", names(colData(object)))
 
})
