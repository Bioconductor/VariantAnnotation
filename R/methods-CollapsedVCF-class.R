### =========================================================================
### CollapsedVCF class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###
### See VCF() on VCF-class page
 
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

### alt() setter
setReplaceMethod("alt", c("CollapsedVCF", "DNAStringSetList"),
    function(x, value)
{
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### expand 
###

setMethod("expand", "CollapsedVCF",
    function(x, ...)
    {
        elt <- elementLengths(alt(x))
        if (all(elt == 1L) || (is(alt(x), "CharacterList"))) {
            fxd <- mcols(fixed(x))[-1]
            fxd$ALT <- unlist(alt(x), use.names=FALSE)
            return(VCF(rowData=rowData(x), colData=colData(x), 
                       exptData=exptData(x), fixed=fxd, 
                       info=mcols(info(x))[-1], geno=geno(x), 
                       ..., collapsed=FALSE))
        }
        idx <- rep.int(seq_len(nrow(x)), elt) 
        hdr <- exptData(x)$header
        ## info
        iexp <- .expandInfo(x, hdr, elt, idx)
        ## fixed
        fexp <- mcols(fixed(x))[idx, -1]
        fexp$ALT <- unlist(alt(x), use.names=FALSE)
        ## geno 
        gexp <- .expandGeno(x, hdr, elt, idx)
        ## rowData
        rdexp <- rowData(x)[idx, ]

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
    icol <- mcols(info(x))[-1]
    inms <- rownames(info(hdr))[info(hdr)$Number == "A"]
    if (length(inms) > 0L) {
        ielt <- sapply(inms, function(i) 
                    elt - elementLengths(mcols(info(x))[[i]]))
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(show, "CollapsedVCF",
    function(object)
{
    .showVCFSubclass(object)
})
