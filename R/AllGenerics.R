setGeneric("predictCoding", signature = c("query", "subject"),
           function(query, subject, seqSource, refAllele, varAllele, ...)
           standardGeneric("predictCoding")
)

setGeneric("getTranscriptSeqs", signature = c("query", "subject"),
           function(query, subject, ...)
           standardGeneric("getTranscriptSeqs")
)

setGeneric("name", function(x, ...) standardGeneric("name"))

setGeneric("stats", function(x, ...) standardGeneric("stats"))

setGeneric("subset", function(x, ...) standardGeneric("subset"))

setGeneric("vaFilter", function(fun, name=NA_character_, ...)
           standardGeneric("vaFilter"),
           signature="fun")


