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
    if (length(gvar) == 0L)
        return(gvar)
    ghdr <- geno(hdr)[rownames(geno(hdr)) %in% names(gvar),]
    ## expand length of ALT
    isA <- ghdr$Number == "A"
    if (any(isA)) {
        gnms <- rownames(ghdr)[ghdr$Number == "A"]
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
    ## list length of ALT each with one REF,ALT pair
    isAD <- names(gvar) == "AD"
    if (any(isAD))
        gvar$AD <- .expandAD(gvar$AD, length(idx), ncol(x)) 
    gvar[!isA & !isAD] <- endoapply(gvar[!isA & !isAD], function(i) {
                              if (is(i, "matrix")) {
                                  matrix(i[idx, ], ncol=ncol(x))
                              } else {
                                  idim <- c(length(idx), dim(i)[2], dim(i)[3])
                                  array(i[idx, , ], dim=idim)
                              }})
    gvar
}

.expandAD <- function(AD, idxlen, xcols)
{
    adpart <- PartitioningByWidth(AD)
    nalt <- width(adpart) - 1L
    if (sum(nalt) != idxlen*xcols)
      stop("length of AD does not match expanded index")
    AD <- unlist(AD, use.names=FALSE)
    ref <- logical(length(AD))
    ref[start(adpart)] <- TRUE
    ##ref <- c(TRUE, tail(seq_along(AD), -1L) %in% (cumsum(elen) + 1L))
    vec <- c(rep(AD[ref], nalt), AD[!ref])
    array(vec, c(idxlen, xcols, 2L))
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

setMethod("expand", "ExpandedVCF",
    function(x, ...)
    {
        x
    }
)

