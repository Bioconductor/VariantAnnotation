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
    if (is.null(names(cdsbytx)) || is.null(names(exbytx)))
        stop("cdsbytx and exbytx objects must have names on the ",
             "outer list elements (transcript names)") 

    ## FIXME : .setActiveSubjectSeq equivalent for GRangesLists
    cdsGR <- unlist(cdsbytx, use.names=FALSE)
    cdsFO <- findOverlaps(ranges, cdsGR, type="within")
    cdsid <- values(unlist(cdsbytx, use.names=FALSE))[["cds_id"]][subjectHits(cdsFO)]
    if (length(cdsFO) == 0)
        return(GRanges())
    nstrand <- as.vector(strand(cdsGR)[subjectHits(cdsFO)] == "-")
    qsub <- ranges[queryHits(cdsFO)]
    txid <- rep(names(cdsbytx), elementLengths(cdsbytx))[subjectHits(cdsFO)]
    txsub <- exbytx[match(txid, names(exbytx))]
    exGR <- unlist(txsub, use.names=FALSE)
    exFO <- findOverlaps(ranges, exGR, type="within")
    c1 <- data.frame(qsub=queryHits(cdsFO), txid=txid)
    tt <- rep(txid, elementLengths(txsub))
    yy <- seq_len(length(exGR))[subjectHits(exFO)]
    c2 <- data.frame(qsub=queryHits(exFO), txid=tt[subjectHits(exFO)],
        exGRindex=yy)
    dups <- duplicated(c2[,c("qsub", "txid")])
    c2 <- c2[!dups,]
    p1 <- paste(c1$qsub, c1$txid, sep=",")
    p2 <- paste(c2$qsub, c2$txid, sep=",")
    c2 <- c2[match(p1, p2), ] 

    cdna <- .refLocsTocDNALocs(qsub, nstrand, txsub, exGR, c2$exGRindex)
    cds <- .refLocsToCDSLocs(qsub, nstrand, cdsbytx, cdsGR, cdsFO)
    p <- c(ceiling(start(cds)/3), ceiling(end(cds)/3))
    protein <- unique(IntegerList(split(p, rep(seq_len(length(p)/2)), 2)))

    values(qsub) <- append(values(qsub), DataFrame(cDNA_loc=cdna, cds_loc=cds, 
        protein_loc=protein, tx_id=txid, cds_id=cdsid))
    qsub
})

.refLocsTocDNALocs <- function(reflocs, nstrand, grlist, lform, olaps)
{
    bounds <- ranges(lform)[olaps]
    cumsums <- .listCumsumShifted(width(grlist))
    qrngs <- ranges(reflocs)
    if (any(nstrand == FALSE))
        qrngs[!nstrand] <- shift(qrngs[!nstrand], - start(bounds)[!nstrand])
    if (any(nstrand))
        qrngs[nstrand] <- IRanges(end(bounds)[nstrand] - end(qrngs)[nstrand],
            width=width(qrngs)[nstrand])

    shift(qrngs, 1L + cumsums[olaps])
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

