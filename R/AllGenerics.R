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
    function(file, genome, param, ...)
    standardGeneric("readVcf")
)

setGeneric("writeVcf", signature = c("obj", "filename"),
    function(obj, filename, ...)
    standardGeneric("writeVcf")
)

## scanVcf

setGeneric("ScanVcfParam",
           function(info=character(), geno=character(), trimEmpty=TRUE,
                    which, asGRanges=character(), ...)
           standardGeneric("ScanVcfParam"),
           signature="which")

setGeneric("scanVcfHeader",
           function(file, ...) standardGeneric("scanVcfHeader"))

setGeneric("scanVcf",
           function(file, ..., param) standardGeneric("scanVcf"))

## VAFilter and VAFilterResult 

setGeneric("name", 
    function(x, ...) 
    standardGeneric("name")
)

setGeneric("stats", 
    function(x, ...) 
    standardGeneric("stats")
)

setGeneric("vaFilter", signature = "fun",
    function(fun, name=NA_character_, ...)
    standardGeneric("vaFilter")
)

setGeneric(".vaValidity")

## VCF class

setGeneric("info", signature = c("x", "i"), 
    function(x, i, ..., withDimnames=TRUE) 
    standardGeneric("info")
)

setGeneric("info<-", signature = c("x", "i", "value"),
    function(x, i, ..., value) 
    standardGeneric("info<-")
)

setGeneric("geno", signature = c("x", "i"),
    function(x, i, ..., withDimnames=TRUE) 
    standardGeneric("geno"),
)

setGeneric("geno<-", signature = c("x", "i", "value"),
    function(x, i, ..., value) 
    standardGeneric("geno<-")
)
















