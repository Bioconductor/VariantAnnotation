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
    isA <- geno(hdr)$Number == "A"
    geno <- geno(x)
    if (any(isA)) {
        gnms <- rownames(geno(hdr))[geno(hdr)$Number == "A"]
        gelt <- sapply(gnms, function(i) 
                    elt - elementLengths(geno[[i]]))
        ## elementLengths same as ALT
        csums <- colSums(gelt) == 0L
        if (any(csums))
            geno[gnms[csums]] <- endoapply(geno[gnms[csums]], function(i)
                                     matrix(unlist(i, use.names=FALSE)))
        ## elementLengths shorter than ALT
        if (any(!csums)) {
            nms <- names(geno) %in% names(csums)[!csums]
            reps <- lapply(list(gelt[!csums] + 1L), rep.int,
                        x=seq_len(nrow(x)))
            geno[nms] <- mendoapply(geno[nms], function(d, r)
                             unlist(d[r], use.names=FALSE),
                         r=reps)
        }
        geno
    }
    geno[!isA] <- endoapply(geno[!isA], function(i) {
                      if (is(i, "matrix"))
                          matrix(i[idx, ], ncol=ncol(x))
                      else
                          i[idx, , ]})
    geno
}

.expandInfo <- function(x, hdr, elt, idx)
{
    icol <- info(x)
    inms <- rownames(info(hdr))[info(hdr)$Number == "A"]
    if (length(inms) > 0L) {
        ielt <- sapply(inms, function(i) 
                    elt - elementLengths(info(x)[[i]]))
        ## elementLengths same as ALT
        csums <- colSums(ielt) == 0L
        if (any(csums))
            res <- expand(icol, inms[csums], TRUE)
        else
            res <- icol[idx, ] 
        ## elementLengths shorter than ALT
        if (any(!csums)) {
            nms <- colnames(icol) %in% names(csums)[!csums]
            reps <- lapply(list(ielt[!csums] + 1L), rep.int,
                        x=seq_len(nrow(x)))
            res[nms] <- Map(function(d, r)
                             unlist(d[r], use.names=FALSE),
                         icol[nms], reps)
        }
        res
    } else {
        icol[idx, ]
    }
}

setMethod("expand", "ExpandedVCF",
    function(x, ...)
    {
        x
    }
)

