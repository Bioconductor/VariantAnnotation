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
        if (all(elt == 1L) || is(alt(x), "CharacterList")) {
            fxd <- mcols(fixed(x))[-1]
            fxd$ALT <- unlist(alt(x), use.names=FALSE)
            return(VCF(rowData=rowData(x), colData=colData(x), 
                       exptData=exptData(x), fixed=fxd, 
                       info=mcols(info(x))[-1], geno=geno(x), 
                       ..., collapsed=FALSE))
        }
        idx <- rep.int(seq_len(nrow(x)), elt) 

        ## info and fixed:
        iexp <- .expandInfo(x, elt, idx)
        fexp <- mcols(fixed(x))[idx, -1]
        fexp$ALT <- unlist(alt(x), use.names=FALSE)

        ## geno: 
        gexp <- endoapply(geno(x), function(i) {
                    if (is(i, "matrix"))
                        i[idx, ]
                    else
                        i[idx, , ]
                 })

        ## rowData:
        rdexp <- rowData(x)[idx, ]

        ## exptData, colData untouched
        VCF(rowData=rdexp, colData=colData(x), exptData=exptData(x),
            fixed=fexp, info=iexp, geno=gexp, ..., collapsed=FALSE)
    }
)

.expandInfo <- function(x, elt, idx)
{
    hd <- exptData(x)$header
    icol <- mcols(info(x))[-1]
    ## FIXME: handle inconsistent elementLengths in readVcf()? 
    inms <- rownames(info(hd))[info(hd)$Number == "A"]
    ielt <- sapply(inms, function(i) 
                elt - elementLengths(mcols(info(x))[[i]]))
    if (any(colSums(ielt) == 0L))
        res <- expand(icol, inms[colSums(ielt) == 0L], TRUE)
    if (any(colSums(ielt) != 0)) {
        res <- icol[idx, ]
        for (i in inms[colSums(ielt) != 0L]) {
            irep <- rep(seq_len(nrow(icol)), ielt[colnames(ielt) == i] + 1L)
            res[[i]] <- unlist(icol[[i]][irep], use.names=FALSE)
        }
    }
    res
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(show, "CollapsedVCF",
    function(object)
{
    .showVCFSubclass(object)
})
