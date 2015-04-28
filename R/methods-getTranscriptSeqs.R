### =========================================================================
### getTranscriptSeqs methods 
### =========================================================================

setMethod("getTranscriptSeqs",  c("GRangesList", "ANY"),
    function(query, subject, ...)
        .Defunct("extractTranscriptSeqs")
)

setMethod("getTranscriptSeqs",  c("GRangesList", "FaFile"),
    function(query, subject, ...)
        .Defunct("extractTranscriptSeqs")
)
 
setMethod("getTranscriptSeqs",  c("GRanges", "FaFile"),
    function(query, subject, ...)
        .Defunct("extractTranscriptSeqs")
)

