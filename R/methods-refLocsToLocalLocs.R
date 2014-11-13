### =========================================================================
### refLocsToLocalLocs methods 
### =========================================================================

setMethod("refLocsToLocalLocs", 
    signature("GRanges", "TxDb", "missing"),
    function(ranges, txdb, cdsbytx, ...)
{
    .Defunct(msg=paste0("refLocsToLocalLocs is deprecated. See ",
                        "?mapCoords methods in GenomicRanges and ",
                        "GenomicAlignments packages.")) 
})

setMethod("refLocsToLocalLocs", 
    signature("GRanges", "missing", "GRangesList"),
    function(ranges, txdb, cdsbytx, ...)
{
    .Defunct(msg=paste0("refLocsToLocalLocs is deprecated. See ",
                        "?mapCoords methods in GenomicRanges and ",
                        "GenomicAlignments packages.")) 
})

