setGeneric("predictCoding", signature = c("query", "subject"),
           function(query, subject, variantColumn, BSgenomeOrganism, ...)
           standardGeneric("predictCoding")
)

setGeneric("getTranscriptSeqs", signature = c("seqSource", "subject"),
           function(seqSource, subject, ...)
           standardGeneric("getTranscriptSeqs")
)
