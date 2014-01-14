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
    ## cds and protein
    cdsGR <- unlist(cdsbytx, use.names=FALSE)
    cdsFO <- findOverlaps(ranges, cdsGR, type="within")
    cdsid <- 
        values(unlist(cdsbytx, use.names=FALSE))[["cds_id"]][subjectHits(cdsFO)]
    if (is.null(cdsid))
        cdsid <- NA_character_
    if (length(cdsFO) == 0)
        return(GRanges())
    nstrand <- as.vector(strand(cdsGR)[subjectHits(cdsFO)] == "-")
    qsub <- ranges[queryHits(cdsFO)]
    names(qsub) <- names(ranges)[queryHits(cdsFO)]
    cds <- .refLocsToCDSLocs(qsub, nstrand, cdsbytx, cdsGR, cdsFO)
    pends <- c(ceiling(start(cds)/3), ceiling(end(cds)/3))
    protein <- unique(IntegerList(split(pends, rep(seq_len(length(pends)/2)), 2)))

    ## set strand to match subject overlap
    strand(qsub) <- strand(cdsGR[subjectHits(cdsFO)]) 
    txid <- rep(names(cdsbytx), elementLengths(cdsbytx))[subjectHits(cdsFO)]
    if (is.null(txid))
        txid <- NA_integer_
    values(qsub) <- append(values(qsub), 
        DataFrame(CDSLOC=cds, 
                  PROTEINLOC=protein, 
                  QUERYID=queryHits(cdsFO), 
                  TXID=txid, CDSID=cdsid))
    qsub
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
