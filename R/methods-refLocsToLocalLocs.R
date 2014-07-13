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
    .Deprecated(msg=paste0("refLocsToLocalLocs is deprecated. See ",
                          "?mapCoords methods in GenomicRanges and ",
                          "GenomicAlignments packages.")) 
    .localCoordinates(ranges, cdsbytx, ...)
})

.localCoordinates <- function(from, to, ...)
{
    ## cds-centric 
    map <- mapCoords(from, to, eltHits=TRUE, ...)
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
        values(unlist(to, use.names=FALSE))[["cds_id"]][eolap]
    if (is.null(cdsid))
        cdsid <- NA_character_
    txid <- rep(names(to), elementLengths(to))[eolap]
    if (is.null(txid))
        txid <- NA_integer_

    ## protein-centric
    res <- from[qolap]
    pends <- c(ceiling(start(map)/3), ceiling(end(map)/3))
    plocs <- unique(IntegerList(split(pends, rep(seq_len(length(pends)/2)), 2)))

    mcols(res) <- append(values(res), 
        DataFrame(CDSLOC=ranges(map), 
                  PROTEINLOC=plocs, 
                  QUERYID=qolap, 
                  TXID=txid, CDSID=cdsid))
    res 
}
