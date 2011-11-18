setGeneric("predictCoding", 
           signature=c("query", "subject", "seqSource", "varAllele"),
           function(query, subject, seqSource, varAllele, ...)
           standardGeneric("predictCoding")
)

setGeneric("locateVariants", signature = c("query", "subject"),
           function(query, subject, ...)
           standardGeneric("locateVariants")
)

setGeneric("getTranscriptSeqs", signature = c("query", "subject"),
           function(query, subject, ...)
           standardGeneric("getTranscriptSeqs")
)

setGeneric("name", function(x, ...) standardGeneric("name"))

setGeneric("stats", function(x, ...) standardGeneric("stats"))

setGeneric("vaFilter", signature = "fun",
           function(fun, name=NA_character_, ...)
           standardGeneric("vaFilter")
)

setGeneric("readVcf", signature = c("file", "genome", "param"),
           function(file, genome, ..., param)
           standardGeneric("readVcf")
)
