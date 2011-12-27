### =========================================================================
### VCF class methods 
### =========================================================================


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Constructor 
##

VCF <-
    function(geno=SimpleList(), rowData=GRanges(), colData=DataFrame(), 
             exptData=SimpleList(), info=SimpleList(),
             ..., verbose=FALSE)
{
    ## FIXME: strip info (row)names
    sx <- SummarizedExperiment(assays=geno, rowData=rowData,
                               colData=colData, exptData=exptData,
                               verbose=verbose)
    new("VCF", sx, info=info, ...)
 }
 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Getters and Setters
##

## info
setMethod(info, c("VCF", "missing"),
    function(x, i, ..., withDimnames=TRUE)
{
    if (withDimnames) {
        endoapply(slot(x, "info"), function(elt, nms) {
            if (is.null(dim(elt))) {
                names(elt) <- nms
                elt
            } else {
                rownames(elt) <- nms
                elt
           }
        }, nms=rownames(x))
    } 
    else
        slot(x, "info")
})

setMethod(info, c("VCF", "numeric"),
    function(x, i, ..., withDimnames=TRUE)
{
    msg <- 'info(<VCF>, i="numeric", ...) invalid subscript "i"'
    tryCatch({
        info(x, ...)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
})

setMethod(info, c("VCF", "character"),
    function(x, i = names(x)[1], ..., withDimnames=TRUE)
{
    msg <- 'info(<VCF>, i="character", ...) invalid subscript "i"'
    res <-
        tryCatch({
            info(x, ...)[[i]]
        }, error=function(err) {
            stop(msg, "\n", conditionMessage(err))
        })
    if (is.null(res))
        stop(msg, "\n    '", i, "' not in names(info(<VCF>))")
    res
})

setReplaceMethod("info", c("VCF", "missing", "SimpleList"),
    function(x, i, ..., value)
{
    initialize(x, ..., info=value)
})

setReplaceMethod("info", c("VCF", "missing", "list"),
    function(x, i, ..., value)
{
    value <- lapply(x, "names<-", NULL)
    initialize(x, ..., info=SimpleList(value))
})

setReplaceMethod("info", c("VCF", "numeric", "Vector"),
    function(x, i = 1, ..., value)
{
    info(x, ...)[[i]] <- value
    x
})

setReplaceMethod("info", c("VCF", "character", "Vector"),
    function(x, i = names(x)[1], ..., value)
{
    info(x, ...)[[i]] <- value
    x
})


## geno
setMethod(geno, c("VCF", "missing"),
    function(x, i, ..., withDimnames=TRUE)
{
    assays(x, ..., withDimnames=withDimnames) 
})

setMethod(geno, c("VCF", "ANY"),
    function(x, i, ..., withDimnames=TRUE)
{
    assay(x, i, ...) 
})

setReplaceMethod("geno", c("VCF", "missing", "ANY"),
    function(x, i, ..., value)
{
    assays(x, ..., value=value)
})

setReplaceMethod("geno", c("VCF", "ANY", "ANY"),
    function(x, i, ..., value)
{
    assays(x, i, ..., value=value)
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

    assays <- endoapply(assays(x, withDimnames=FALSE), function(elt) {
                  if (class(elt) == "array")
                     elt[i, j, , drop=FALSE]
                  else
                     elt[i, j, drop=FALSE]
              })
    info <- endoapply(info(x), function(elt) {
                  if (class(elt) == "array")
                     elt[i, , , drop=FALSE]
                  else
                     elt[i, , drop=FALSE]
              })
 
    initialize(x, 
               rowData=rowData(x)[i,,drop=FALSE],
               colData=colData(x)[j,,drop=FALSE],
               assays=assays, info=info)
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
                   nms <- names(r)
                   nms[i] <- names(rowData(value))
                   names(r) <- nms
                   r
               }), colData=local({
                   c <- colData(x)
                   c[j,] <- colData(value)
                   rownames(c)[j] <- rownames(colData(value))
                   c
               }), assays=local({
                   a <- assays(x, withDimnames=FALSE)
                   v <- assays(value, withDimnames=FALSE)
                   mendoapply(function(x, ..., value) {
                      if (class(x) == "array")
                         x[i, j, ] <- value
                      else
                         x[i, j] <- value
                      x
                   }, x=a, value=v, ...)
               }), info=local({
                   ii <- info(x, withDimnames=FALSE)
                   v <- info(value, withDimnames=FALSE)
                   mendoapply(function(x, ..., value) {
                      if (class(x) == "array")
                         x[i, , ] <- value
                      else
                         x[i, ] <- value
                      x
                   }, x=ii, value=v, ...)
               }))
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
    expt <- names(exptData(object))
    if (is.null(expt))
        expt <- character(length(exptData(object)))
    scat("exptData(%d): %s\n", expt)
    info <- names(info(object, withDimnames=FALSE))
    if (is.null(info))
        info <- character(length(info(object, withDimnames=FALSE)))
    scat("info(%d): %s\n", info)
    nms <- names(assays(object, withDimnames=FALSE))
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
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
