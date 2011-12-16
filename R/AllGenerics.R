setGeneric("predictCoding", 
           signature=c("query", "subject", "seqSource", "varAllele"),
           function(query, subject, seqSource, varAllele, ...)
           standardGeneric("predictCoding")
)

setGeneric("locateVariants", signature = c("query", "subject"),
           function(query, subject, ...)
           standardGeneric("locateVariants")
)

setGeneric("MatrixToSnpMatrix", signature = c("callMatrix", "mapGRanges"),
           function(callMatrix, mapGRanges, ...)
           standardGeneric("MatrixToSnpMatrix")
)

setGeneric("getTranscriptSeqs", signature = c("query", "subject"),
           function(query, subject, ...)
           standardGeneric("getTranscriptSeqs")
)

setGeneric("readVcf", signature = c("file", "genome", "param"),
           function(file, genome, ..., param)
           standardGeneric("readVcf")
)

setGeneric("writeVcf", signature = c("file", "obj"),
           function(file, obj, ...)
           standardGeneric("writeVcf")
)

setGeneric("name", function(x, ...) standardGeneric("name"))

setGeneric("stats", function(x, ...) standardGeneric("stats"))

setGeneric("vaFilter", signature = "fun",
           function(fun, name=NA_character_, ...)
           standardGeneric("vaFilter")
)

