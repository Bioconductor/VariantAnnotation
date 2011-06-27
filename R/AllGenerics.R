setGeneric("predictCoding", signature = c("query", "subject"),
           function(query, subject, seqSource, refAllele, varAllele, ...)
           standardGeneric("predictCoding")
)

setGeneric("getTranscriptSeqs", signature = c("query", "subject"),
           function(query, subject, ...)
           standardGeneric("getTranscriptSeqs")
)


### SRFilter
#
#setGeneric("name", function(x, ...) standardGeneric("name"))
#
#setGeneric("stats", function(x, ...) standardGeneric("stats"))
#
#setGeneric("srFilter", function(fun, name=NA_character_, ...)
#           standardGeneric("srFilter"),
#           signature="fun")
#

