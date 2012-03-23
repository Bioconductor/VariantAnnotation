### =========================================================================
### refLocsToLocalLocs methods 
### =========================================================================

setMethod("refLocsToLocalLocs", 
    signature("GRanges", "TranscriptDb", "missing", "missing"),
    function(ranges, txdb, cdsbytx, exbytx, ...)
{
    ## remove circular
    .setActiveSubjectSeq(ranges, txdb)
    cdsbytx <- cdsBy(txdb, "tx")
    exbytx <- exonsBy(txdb, "tx")
    callGeneric(ranges=ranges, cdsbytx=cdsbytx, exbytx=exbytx)
})

setMethod("refLocsToLocalLocs", 
    signature("GRanges", "missing", "GRangesList", "GRangesList"),
    function(ranges, txdb, cdsbytx, exbytx, ...)
{
    ## FIXME : .setActiveSubjectSeq equivalent for GRangesLists
    if (is.null(names(cdsbytx)) || is.null(names(exbytx)))
        stop("the outer list elements of cdsbytx and exbytx must have ",
             " names (i.e., transcript identifiers)") 

    ## cds and protein
    cdsGR <- unlist(cdsbytx, use.names=FALSE)
    cdsFO <- findOverlaps(ranges, cdsGR, type="within")
    cdsid <- values(unlist(cdsbytx, use.names=FALSE))[["cds_id"]][subjectHits(cdsFO)]
    if (length(cdsFO) == 0)
        return(GRanges())
    nstrand <- as.vector(strand(cdsGR)[subjectHits(cdsFO)] == "-")
    qsub <- ranges[queryHits(cdsFO)]
    cds <- .refLocsToCDSLocs(qsub, nstrand, cdsbytx, cdsGR, cdsFO)
    pends <- c(ceiling(start(cds)/3), ceiling(end(cds)/3))
    protein <- unique(IntegerList(split(pends, rep(seq_len(length(pends)/2)), 2)))

    ## cDNA
    txid <- rep(names(cdsbytx), elementLengths(cdsbytx))[subjectHits(cdsFO)]
    exbytx <- exbytx[match(txid, names(exbytx))]
    exGR <- unlist(exbytx, use.names=FALSE)
    exFO <- findOverlaps(ranges, exGR, type="within")
    cdna <- .refLocsTocDNALocs(qsub, nstrand, exbytx, exGR, exFO, cdsFO)

    values(qsub) <- append(values(qsub), DataFrame(cDNA_loc=cdna, cds_loc=cds, 
        protein_loc=protein, tx_id=txid, cds_id=cdsid))
    qsub
})

.refLocsTocDNALocs <- function(reflocs, nstrand, grlist, lform, olaps, cdsolaps)
{
    ## keep one record per query-tx hit
    tx <- rep(names(grlist), elementLengths(grlist))[subjectHits(olaps)]
    idx <- seq_len(length(lform))[subjectHits(olaps)]
    map <- data.frame(query=queryHits(olaps), tx=tx, idx=idx)
    dups <- duplicated(map[,c("query", "tx")])
    map <- map[!dups,]

    ## remove ranges in exons outside cds regions
    cds <- map$query %in% unique(queryHits(cdsolaps))
    map <- map[cds,]

    bounds <- ranges(lform)[map$idx]
    cumsums <- .listCumsumShifted(width(grlist))
    qrngs <- ranges(reflocs)
    if (any(nstrand == FALSE))
        qrngs[!nstrand] <- shift(qrngs[!nstrand], - start(bounds)[!nstrand])
    if (any(nstrand))
        qrngs[nstrand] <- IRanges(end(bounds)[nstrand] - end(qrngs)[nstrand],
            width=width(qrngs)[nstrand])

    shift(qrngs, 1L + cumsums[map$idx])
}

.refLocsToCDSLocs <- function(reflocs, nstrand, grlist, lform, olaps)
{
    bounds <- ranges(lform)[subjectHits(olaps)]
    cumsums <- .listCumsumShifted(width(grlist))
    qrngs <- ranges(reflocs)
    if (any(nstrand == FALSE))
        qrngs[!nstrand] <- shift(qrngs[!nstrand], - start(bounds)[!nstrand])
    if (any(nstrand))
        qrngs[nstrand] <- IRanges(end(bounds)[nstrand] - end(qrngs)[nstrand],
            width=width(qrngs)[nstrand])

    shift(qrngs, 1L + cumsums[subjectHits(olaps)])
}

