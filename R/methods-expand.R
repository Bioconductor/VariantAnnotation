### =========================================================================
### expand methods 
### =========================================================================

setMethod("expand", "CollapsedVCF",
    function(x, ..., row.names = FALSE)
    {
        rd <- rowData(x, fixed=FALSE)
        if (!row.names)
            names(rd) <- NULL

        elt <- elementLengths(alt(x))
        if (all(elt == 1L)) {
            fxd <- fixed(x)
            fxd$ALT <- unlist(alt(x), use.names=FALSE)
            AD <- "AD" %in% names(geno(x))
            if (AD)
              geno(x)$AD <- .expandAD(geno(x)$AD, seq_len(nrow(x)), ncol(x))
            return(VCF(rd, colData(x), exptData(x), fxd, .unlistAltInfo(x), 
                       geno(x), ..., collapsed=FALSE))
        }

        ## info, fixed, geno
        hdr <- exptData(x)$header
        idx <- rep.int(seq_len(nrow(x)), elt)
        iexp <- .expandInfo(x, hdr, elt, idx)
        fexp <- fixed(x)[idx, ]
        fexp$ALT <- unlist(alt(x), use.names=FALSE)
        gexp <- .expandGeno(x, hdr, elt, idx)

        ## rowData
        if (is.null(rd$paramRangeID)) {
            rdexp <- rd[idx, ]
            mcols(rdexp) <- NULL
        } else rdexp <- rd[idx, "paramRangeID"]

        VCF(rdexp, colData(x), exptData(x), fexp, iexp, gexp,
            ..., collapsed=FALSE)
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
        gelt <- do.call(cbind, lapply(gnms, function(i) 
                                      elt - elementLengths(gvar[[i]])))
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
                if ((length(unique(elt)) != 1L) && 
                    (dim(AD)[3] != unique(elt))) {
                    warning("'AD' was ignored: number of 'AD' values ",
                            "do not match REF + ALT")
                    isAD[isAD] <- FALSE 
                }
            } else {
               isAD[isAD] <- FALSE 
            }
        } else {
            ## 'Number' is '.'
            gvar$AD <- .expandAD(AD, idx, ncol(x))
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

## returns an array of REF,ALT,<NON_REF> pairs
.expandAD <- function(AD, idx, xcols)
{
    if (is.list(AD)) {
        part <- PartitioningByWidth(AD)
        maxelt <- max(width(part))
        if (any(zeros <- width(part) == 0L)) 
            AD[zeros] <- list(rep(NA_integer_, maxelt))
        AD <- AD[idx]
        array(t(do.call(cbind, AD)), c(length(idx), xcols, maxelt))
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
    inms <- names(ivar)[info(hdr)[names(ivar), "Number"] == "A"]
    if (length(inms) > 0L) {
        ielt <- do.call(cbind, lapply(inms, function(i) 
                                      elt - elementLengths(ivar[[i]])))
        ## elementLengths same as ALT
        csums <- colSums(ielt) == 0L
        if (any(csums))
            res <- IRanges:::.expandByColumnSet(ivar, inms[csums], TRUE)
        else
            res <- ivar[idx, , drop=FALSE] 
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
        ivar[idx, , drop=FALSE]
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
    function(x, ...) x
)

