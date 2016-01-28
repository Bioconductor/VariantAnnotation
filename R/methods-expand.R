### =========================================================================
### expand methods 
### =========================================================================

setMethod("expand", "CollapsedVCF",
    function(x, ..., row.names = FALSE)
    {
        rd <- rowRanges(x, fixed=FALSE)
        if (!row.names)
            names(rd) <- NULL

        elt <- elementLengths(alt(x))
        if (all(elt == 1L)) {
            fxd <- fixed(x)
            fxd$ALT <- unlist(alt(x), use.names=FALSE)
            AD <- "AD" %in% names(geno(x))
            if (AD)
              geno(x)$AD <- .expandAD(geno(x)$AD, nrow(x), ncol(x))
            return(VCF(rd, colData(x), metadata(x), fxd, .unlistAltInfo(x), 
                       geno(x), ..., collapsed=FALSE))
        }

        ## info, fixed, geno
        hdr <- metadata(x)$header
        idx <- rep.int(seq_len(nrow(x)), elt)
        iexp <- .expandInfo(x, hdr, elt, idx)
        fexp <- fixed(x)[idx, ]
        fexp$ALT <- unlist(alt(x), use.names=FALSE)
        gexp <- .expandGeno(x, hdr, elt, idx)

        ## rowRanges
        if (is.null(rd$paramRangeID)) {
            rdexp <- rd[idx, ]
            mcols(rdexp) <- NULL
        } else rdexp <- rd[idx, "paramRangeID"]

        VCF(rdexp, colData(x), metadata(x), fexp, iexp, gexp,
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

## returns an array of REF,ALT and REF,<NON_REF> pairs
## 'z' dimension of result is always 2
.expandAD <- function(AD, idxlen, xcols)
{
    if (is.list(AD)) {
        adpart <- PartitioningByWidth(AD)
        if (any(zeros <- width(adpart) == 0L)) { 
            AD[zeros] <- list(rep(NA_integer_, 2L))
            adpart <- PartitioningByWidth(AD)
        }

        nalt <- width(adpart) - 1L
        AD <- unlist(AD, use.names=FALSE)
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
    hnms <- intersect(names(ivar), rownames(info(hdr)))
    inms <- hnms[info(hdr)[hnms, "Number"] == "A"]
    if (length(inms) > 0L) {
        ielt <- do.call(cbind, lapply(inms, function(i) 
                                      elt - elementLengths(ivar[[i]])))
        ## elementLengths same as ALT
        csums <- colSums(ielt) == 0L
        if (any(csums))
            res <- S4Vectors:::expandByColumnSet(ivar, inms[csums], TRUE)
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

