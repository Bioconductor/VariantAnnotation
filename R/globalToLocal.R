globalToLocal <- function(query, subject) 
{
    ## convert genome coordinates to transcript coordinates
    ## complement of transcriptLocs2refLocs() in GenomicFeatures 
    ## subject must be a cdsbytx object
    gr <- unlist(subject, use.names=FALSE)
    fo <- findOverlaps(query, gr, type = "within")
    shit <- subjectHits(fo)
    qhit <- queryHits(fo)

    ## position of query in coding region 
    qrng <- ranges(query)[qhit]
    cds_bounds <- ranges(gr)[shit]
    neg <- as.vector(strand(gr)[shit] == "-")
    if (any(neg == FALSE))
        qrng[!neg] <- shift(qrng[!neg], - start(cds_bounds)[!neg])
    if (any(neg))
        qrng[neg] <- IRanges(end(cds_bounds)[neg] - end(qrng)[neg],
            width=width(qrng)[neg])
    cdsloc <- start(qrng)
    strand <- rep("+", length(neg))
    strand[neg] <- "-"

    ## position of query in transcript 
    cumsums <- .listCumsumShifted(width(subject))
    txloc <- shift(qrng, 1L + cumsums[shit])

    sindex <- rep(seq(length(subject)), elementLengths(subject))[shit]
    cdsid <- values(gr[shit])[["cds_id"]]
    DataFrame(qindex=qhit, sindex, strand, cdsid, txloc, cdsloc)
}

