setGeneric("predictCoding", signature = c("query", "subject"),
           function(query, subject, variantAllele, seqSource, ...)
           standardGeneric("predictCoding")
)

setGeneric("getTranscriptSeqs", signature = c("query", "subject"),
           function(query, subject, ...)
           standardGeneric("getTranscriptSeqs")
)
