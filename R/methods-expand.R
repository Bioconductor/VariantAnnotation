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
            if (!is.null(geno(x)$AD))
              geno(x)$AD <- .expandAD(geno(x)$AD, nrow(x), ncol(x))
            return(VCF(rowData=rowData(x), colData=colData(x), 
                       exptData=exptData(x), fixed=fxd, 
                       info=.unlistAltInfo(x), geno=geno(x), 
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
        if (is.null(rowData(x)$paramRangeID)) {
            rdexp <- rowData(x)[idx, ]
            mcols(rdexp) <- NULL
        } else {
            rdexp <- rowData(x)[idx, "paramRangeID"]
        }
        ## exptData, colData untouched
        VCF(rowData=rdexp, colData=colData(x), exptData=exptData(x),
            fixed=fexp, info=iexp, geno=gexp, ..., collapsed=FALSE)
    }
)

.expandGeno <- function(x, hdr, elt, idx)
{
    gvar <- geno(x)
    if (length(gvar) == 0L)
        return(gvar)
    ghdr <- geno(hdr)[rownames(geno(hdr)) %in% names(gvar),]
    ## 'Number=A': one value per ALT
    isA <- ghdr$Number == "A"
    if (any(isA)) {
        gnms <- rownames(ghdr)[isA]
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
    ## AD field: one value for REF and each ALT
    if (any(isAD <- names(gvar) == "AD")) {
        AD <- gvar$AD
        if (!is.list(AD)) {
            ## 'Number' is integer
            if (is(AD, "array")) {
                if ((length(unique(elt)) != 1L) && (dim(AD)[3] != unique(elt))) {
                    warning("'AD' was ignored: number of 'AD' values ",
                            "do not match REF + ALT")
                    isAD[isAD] <- FALSE 
                }
            } else {
               isAD[isAD] <- FALSE 
            }
        } else {
            ## 'Number' is '.'
            gvar$AD <- .expandAD(AD, length(idx), ncol(x))
        }
    }
    isA <- names(gvar) %in% rownames(ghdr)[isA]
    gvar[!isA & !isAD] <- endoapply(gvar[!isA & !isAD], function(i) {
                              if (is(i, "matrix")) {
                                  matrix(i[idx, ], nrow=length(idx),
                                         ncol=ncol(x))
                              } else {
                                  idim <- c(length(idx), dim(i)[2], dim(i)[3])
                                  array(i[idx, , ], dim=idim)
                              }})
    gvar
}

## returns an array of REF,ALT pairs
.expandAD <- function(AD, idxlen, xcols)
{
    if (is.list(AD)) {
        adpart <- PartitioningByWidth(AD)
        nalt <- width(adpart) - 1L
        if (sum(nalt) != idxlen*xcols)
          stop("length of AD does not match expanded index")
        AD <- as.integer(unlist(AD, use.names=FALSE))
        ref <- logical(length(AD))
        ref[start(adpart)] <- TRUE
        vec <- c(rep(AD[ref], nalt), AD[!ref])
        array(vec, c(idxlen, xcols, 2L))
    } else {
        AD
    }
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
                    elt - elementLengths(ivar[[i]]))
        ## elementLengths same as ALT
        csums <- colSums(ielt) == 0L
        if (any(csums))
            res <- IRanges:::.expandByColumnSet(ivar, inms[csums], TRUE)
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

.unlistAltInfo <- function(x) {
    ivar <- info(x)
    if (ncol(ivar) == 0L)
        return(ivar)
    hdr <- header(x)
    inms <- rownames(info(hdr))[info(hdr)$Number == "A"]
    inms <- inms[inms %in% names(ivar)]
    if (length(inms))
        return(ivar)
    ivar[inms] <- lapply(ivar[inms], drop)
    ivar
}

setMethod("expand", "ExpandedVCF",
    function(x, ...)
    {
        x
    }
)

