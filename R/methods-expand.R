### =========================================================================
### expand methods 
### =========================================================================

setMethod("expand", "CollapsedVCF",
    function(x, ...)
    {
        elt <- elementLengths(alt(x))
        if (all(elt == 1L)) {
            fxd <- fixed(x)
            fxd$ALT <- unlist(alt(x), use.names=FALSE)
            return(VCF(rowData=rowData(x), colData=colData(x), 
                       exptData=exptData(x), fixed=fxd, 
                       info=info(x), geno=geno(x), 
                       ..., collapsed=FALSE))
        }
        idx <- rep.int(seq_len(nrow(x)), elt) 
        hdr <- exptData(x)$header
        ## info
        iexp <- .expandInfo(x, hdr, elt, idx)
        ## fixed
        fexp <- fixed(x)[idx, ]
        fexp$ALT <- unlist(alt(x), use.names=FALSE)
        ## geno 
        gexp <- .expandGeno(x, hdr, elt, idx)
        ## rowData
        rdexp <- rowData(x)[idx, "paramRangeID"]

        ## exptData, colData untouched
        VCF(rowData=rdexp, colData=colData(x), exptData=exptData(x),
            fixed=fexp, info=iexp, geno=gexp, ..., collapsed=FALSE)
    }
)

.expandGeno <- function(x, hdr, elt, idx)
{
    gvar <- geno(x)
    ## no data 
    if (length(gvar) == 0L)
        return(gvar)
    isA <- geno(hdr)$Number == "A"
    if (any(isA)) {
        gnms <- rownames(geno(hdr))[geno(hdr)$Number == "A"]
        gelt <- sapply(gnms, function(i) 
                    elt - elementLengths(gvar[[i]]))
        ## elementLengths same as ALT
        csums <- colSums(gelt) == 0L
        if (any(csums))
            gvar[gnms[csums]] <- endoapply(gvar[gnms[csums]], function(i)
                                     matrix(unlist(i, use.names=FALSE)))
        ## elementLengths shorter than ALT
        if (any(!csums)) {
            nms <- names(gvar) %in% names(csums)[!csums]
            reps <- lapply(list(gelt[!csums] + 1L), rep.int,
                        x=seq_len(nrow(x)))
            gvar[nms] <- mendoapply(gvar[nms], function(d, r)
                             unlist(d[r], use.names=FALSE),
                         r=reps)
        }
        gvar
    }
    gvar[!isA] <- endoapply(gvar[!isA], function(i) {
                      if (is(i, "matrix"))
                          matrix(i[idx, ], ncol=ncol(x))
                      else
                          i[idx, , ]})
    gvar
}

.expandInfo <- function(x, hdr, elt, idx)
{
    ivar <- info(x)
    ## no data 
    if (ncol(ivar) == 0L)
        return(DataFrame(row.names=seq_along(idx)))
    inms <- rownames(info(hdr))[info(hdr)$Number == "A"]
    if (length(inms) > 0L) {
        ielt <- sapply(inms, function(i) 
                    elt - elementLengths(info(x)[[i]]))
        ## elementLengths same as ALT
        csums <- colSums(ielt) == 0L
        if (any(csums))
            res <- expand(ivar, inms[csums], TRUE)
        else
            res <- ivar[idx, ] 
        ## elementLengths shorter than ALT
        if (any(!csums)) {
            nms <- colnames(ivar) %in% names(csums)[!csums]
            reps <- lapply(list(ielt[!csums] + 1L), rep.int,
                        x=seq_len(nrow(x)))
            res[nms] <- Map(function(d, r)
                             unlist(d[r], use.names=FALSE),
                         ivar[nms], reps)
        }
        res
    } else {
        ivar[idx, ]
    }
}

setMethod("expand", "ExpandedVCF",
    function(x, ...)
    {
        x
    }
)

