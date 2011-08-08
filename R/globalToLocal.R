globalToLocal <- function(global, ranges) 
{
    gr <- unlist(ranges, use.names=FALSE)
    #gr <- ranges@unlistData
    ol <- findOverlaps(global, gr, type = "within")
    shits <- subjectHits(ol)
    qhits <- queryHits(ol)
    local <- ranges(global)[qhits]
    bounds <- ranges(gr)[shits]
 
    ## location wrt start of coding region 
    neg <- as.vector(strand(gr)[shits] == "-")
    local[!neg] <- shift(local[!neg], - start(bounds)[!neg])
    local[neg] <- IRanges(end(bounds)[neg] - end(local)[neg],
        width = width(local)[neg])

    ## location wrt transcript 
    cumsums <- .cumsumShifted(width(ranges))
    local <- shift(local, 1L + cumsums[shits])

    rangesInd <- rep(seq(length(ranges)), elementLengths(ranges))[shits]
    DataFrame(globalInd = qhits, rangesInd, local)
}

