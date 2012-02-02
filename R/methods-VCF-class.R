### =========================================================================
### VCF class methods 
### =========================================================================


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Constructor 
##

VCF <-
    function(rowData=GRanges(), colData=DataFrame(), exptData=SimpleList(), 
             fixedFields=DataFrame(), info=DataFrame(), geno=SimpleList(),
             ..., verbose=FALSE)
{
    rownames(info) <- rownames(fixedFields) <- NULL
    sx <- SummarizedExperiment(assays=geno, rowData=rowData,
                               colData=colData, exptData=exptData,
                               verbose=verbose)
    new("VCF", sx, info=info, fixedFields=fixedFields, ...)
 }
 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Getters and Setters
##

## ref 
setMethod("ref", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    values(gr) <- DataFrame(REF=slot(x, "fixedFields")$REF)
    gr
})

setReplaceMethod("ref", c("VCF", "DNAStringSet"),
    function(x, value)
{
    slot(x, "fixedFields")$REF <- value
    x
})

## alt 
setMethod("alt", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    values(gr) <- DataFrame(ALT=slot(x, "fixedFields")$ALT)
    gr
})

setReplaceMethod("alt", c("VCF", "CharacterList"),
    function(x, value)
{
    slot(x, "fixedFields")$ALT <- value
    x
})

setReplaceMethod("alt", c("VCF", "DNAStringSetList"),
    function(x, value)
{
    slot(x, "fixedFields")$ALT <- value
    x
})

## qual 
setMethod("qual", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    values(gr) <- DataFrame(QUAL=slot(x, "fixedFields")$QUAL)
    gr
})

setReplaceMethod("qual", c("VCF", "integer"),
    function(x, value)
{
    slot(x, "fixedFields")$QUAL <- value
    x
})

## filt
setMethod("filt", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    values(gr) <- DataFrame(FILTER=slot(x, "fixedFields")$FILTER)
    gr
})

setReplaceMethod("filt", c("VCF", "character"),
    function(x, value)
{
    slot(x, "fixedFields")$FILTER <- value
    x
})

## fixedFields 
setMethod("fixedFields", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    values(gr) <- slot(x, "fixedFields")
    gr
})

setReplaceMethod("filt", c("VCF", "character"),
    function(x, value)
{
    slot(x, "fixedFields") <- value
    x
})

## fixed
setMethod("fixed", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    values(gr) <- DataFrame(slot(x, "fixedFields"), slot(x, "info"))
    gr
})

## info 
setMethod("info", "VCF", 
    function(x) 
{
    gr <- rowData(x)
    values(gr) <- slot(x, "info")
    gr
})

setReplaceMethod("info", c("VCF", "DataFrame"),
    function(x, value)
{
    slot(x, "info") <- value
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
    assay(x, i, ..., value=value)
})

setReplaceMethod("geno", c("VCF", "numeric", "matrix"),
    function(x, i, ..., value)
{
    assay(x, i, ..., value=value)
})

setReplaceMethod("geno", c("VCF", "missing", "SimpleList"),
    function(x, i, ..., value)
{
    assay(x, ..., value=value)
})

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Subsetting 
##

.VCF.subset <-
    function(x, i, j, ...)
{
    if (is.character(i)) {
        msg <- "<VCF>[i,] index out of bounds: %s"
        i <- GenomicRanges:::.SummarizedExperiment.charbound(i, rownames(x), msg)
    }
    if (is.character(j)) {
        msg <- "<VCF>[,j] index out of bounds: %s"
        j <- GenomicRanges:::.SummarizedExperiment.charbound(j, colnames(x), msg)
    }

    geno <- endoapply(geno(x, withDimnames=FALSE), function(elt) {
                  if (class(elt) == "array")
                     elt[i, j, , drop=FALSE]
                  else
                     elt[i, j, drop=FALSE]
              })
 
    initialize(x, 
               rowData=rowData(x)[i,,drop=FALSE],
               colData=colData(x)[j,,drop=FALSE],
               assays=geno, 
               info=values(info(x))[i,,drop=FALSE],
               fixedFields=values(fixedFields(x))[i,])
}

setMethod("[", c("VCF", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,VCF,ANY,ANY-method'")
    if (missing(i) && missing(j))
        x
    else if (missing(i))
        .VCF.subset(x, TRUE, j, ...)
    else if (missing(j))
        .VCF.subset(x, i, TRUE, ...)
    else
        .VCF.subset(x, i, j, ...)
})

.VCF.subsetassign <-
    function(x, i, j, ..., value)
{
    if (is.character(i)) {
        msg <- "<VCF>[i,]<- index out of bounds: %s"
        i <- GenomicRanges:::.SummarizedExperiment.charbound(i, rownames(x), msg)
    }
    if (is.character(j)) {
        msg <- "<VCF>[,j]<- index out of bounds: %s"
        j <- GenomicRanges:::.SummarizedExperiment.charbound(j, colnames(x), msg)
    }
    initialize(x,
               exptData=c(exptData(x), exptData(value)),
               rowData=local({
                   r <- rowData(x)
                   r[i,] <- rowData(value)
                   names(r)[i] <- names(rowData(value))
                   r
               }), colData=local({
                   c <- colData(x)
                   c[j,] <- colData(value)
                   rownames(c)[j] <- rownames(colData(value))
                   c
               }), assays=local({
                   a <- geno(x, withDimnames=FALSE)
                   v <- geno(value, withDimnames=FALSE)
                   mendoapply(function(x, ..., value) {
                      if (class(x) == "array")
                         x[i, j, ] <- value
                      else
                         x[i, j] <- value
                      x
                   }, x=a, value=v, ...)
               }), info=local({
                   ii <- values(info(x))
                   ii[i,] <- values(info(value))
                   ii 
               }), fixedFields=local({
                   ff <- values(fixedFields(x)) 
                   ff[i,] <- values(fixedFields(value)) 
                   ff})
              )
}

setReplaceMethod("[",
    c("VCF", "ANY", "ANY", "VCF"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        x
    else if (missing(i))
        .VCF.subsetassign(x, TRUE, j, ..., value=value)
    else if (missing(j))
        .VCF.subsetassign(x, i, TRUE, ..., value=value)
    else
        .VCF.subsetassign(x, i, j, ..., value=value)
})

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Other methods 
##

setMethod("genome", "VCF",
    function(x)
{
    genome(rowData(x))
})

setMethod("seqlevels", "VCF",
    function(x)
{
    seqlevels(rowData(x))
})

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
    cat("genome:", genome(rowData(object)), "\n")

    expt <- names(exptData(object))
    if (is.null(expt))
        expt <- character(length(exptData(object)))
    scat("exptData(%d): %s\n", expt)

    fixed <- names(values(fixedFields(object)))
    if (is.null(fixed))
        fixed <- character(ncol(values(fixedFields(object))))
    scat("fixedFields(%d): %s\n", fixed)

    info <- names(values(info(object)))
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
