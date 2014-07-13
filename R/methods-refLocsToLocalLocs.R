### =========================================================================
### refLocsToLocalLocs methods 
### =========================================================================

setMethod("refLocsToLocalLocs", 
    signature("GRanges", "TranscriptDb", "missing"),
    function(ranges, txdb, cdsbytx, ...)
{
    cdsbytx <- cdsBy(txdb, "tx")
    callGeneric(ranges=ranges, cdsbytx=cdsbytx)
})

setMethod("refLocsToLocalLocs", 
    signature("GRanges", "missing", "GRangesList"),
    function(ranges, txdb, cdsbytx, ...)
{
    if (is.null(names(cdsbytx)))
        stop("the outer list elements of cdsbytx must have ",
             " names (i.e., transcript identifiers)") 
    ## cds-centric 
    map <- mapCoords(ranges, cdsbytx, eltHits=TRUE, ...)
    if (length(map) == 0) {
        gr <- GRanges()
        mcols(gr) <- DataFrame(REF=DNAStringSet(), ALT=DNAStringSetList(),
                               varAllele=DNAStringSet(), CDSLOC=IRanges(),
                               PROTEINLOC=IntegerList(), QUERYID=integer(),
                               TXID=character(), CDSID=integer())
        return(gr)
    }
    eolap <- map$eltHits
    qolap <- map$queryHits
    cdsid <- 
        values(unlist(cdsbytx, use.names=FALSE))[["cds_id"]][eolap]
    if (is.null(cdsid))
        cdsid <- NA_character_
    txid <- rep(names(cdsbytx), elementLengths(cdsbytx))[eolap]
    if (is.null(txid))
        txid <- NA_integer_

    ## protein-centric
    res <- ranges[qolap]
    pends <- c(ceiling(start(map)/3), ceiling(end(map)/3))
    plocs <- unique(IntegerList(split(pends, rep(seq_len(length(pends)/2)), 2)))

    mcols(res) <- append(values(res), 
        DataFrame(CDSLOC=ranges(map), 
                  PROTEINLOC=plocs, 
                  QUERYID=qolap, 
                  TXID=txid, CDSID=cdsid))
    res 
})

.refLocsToCDSLocs <- function(reflocs, nstrand, grlist, lform, olaps)
{
    bounds <- ranges(lform)[subjectHits(olaps)]
    ## assumption :
    ## cds regions are sorted 5' to 3' (i.e., cds rank is lowest at 5' end)
    cumsums <- .listCumsumShifted(width(grlist))
    qrngs <- ranges(reflocs)
    if (any(nstrand == FALSE))
        qrngs[!nstrand] <- shift(qrngs[!nstrand], - start(bounds)[!nstrand])
    if (any(nstrand))
        qrngs[nstrand] <- IRanges(end(bounds)[nstrand] - end(qrngs)[nstrand],
            width=width(qrngs)[nstrand])

    shift(qrngs, 1L + cumsums[subjectHits(olaps)])
}
