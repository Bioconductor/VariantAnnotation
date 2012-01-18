setGeneric("predictCoding", 
    signature=c("query", "subject", "seqSource", "varAllele"),
    function(query, subject, seqSource, varAllele, ...)
    standardGeneric("predictCoding")
)

setGeneric("locateVariants", signature = c("query", "subject"),
    function(query, subject, ...)
    standardGeneric("locateVariants")
)

setGeneric("MatrixToSnpMatrix", signature = c("callMatrix", "ref", "alt"),
    function(callMatrix, ref, alt, ...)
    standardGeneric("MatrixToSnpMatrix")
)

setGeneric("getTranscriptSeqs", signature = c("query", "subject"),
    function(query, subject, ...)
    standardGeneric("getTranscriptSeqs")
)

## read/write Vcf

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

## VCF class

setGeneric("fixedFields", signature = "x", 
    function(x) standardGeneric("fixedFields")
)

setGeneric("fixedFields<-", signature = c("x", "value"),
    function(x, value) standardGeneric("fixedFields<-")
)

setGeneric("fixed", signature = "x", 
    function(x) standardGeneric("fixed")
)

setGeneric("ref", signature = "x", 
    function(x) standardGeneric("ref")
)

setGeneric("ref<-", signature = c("x", "value"),
    function(x, value) standardGeneric("ref<-")
)

setGeneric("alt", signature = "x", 
    function(x) standardGeneric("alt")
)

setGeneric("alt<-", signature = c("x", "value"),
    function(x, value) standardGeneric("alt<-")
)

setGeneric("qual", signature = "x", 
    function(x) standardGeneric("qual")
)

setGeneric("qual<-", signature = c("x", "value"),
    function(x, value) standardGeneric("qual<-")
)

setGeneric("filt", signature = "x", 
    function(x) standardGeneric("filt")
)

setGeneric("filt<-", signature = c("x", "value"),
    function(x, value) standardGeneric("filt<-")
)

setGeneric("info", signature = "x", 
    function(x) standardGeneric("info")
)

setGeneric("info<-", signature = c("x", "value"),
    function(x, value) standardGeneric("info<-")
)

setGeneric("geno", signature = c("x", "i"),
    function(x, i, ..., withDimnames=TRUE) 
    standardGeneric("geno"),
)

setGeneric("geno<-", signature = c("x", "i", "value"),
    function(x, i, ..., value) 
    standardGeneric("geno<-")
)

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
















