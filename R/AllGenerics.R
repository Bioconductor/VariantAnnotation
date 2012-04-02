### -------------------------------------------------------------------------
### predictCoding 
###

setGeneric("predictCoding", 
    signature=c("query", "subject", "seqSource", "varAllele"),
    function(query, subject, seqSource, varAllele, ...)
    standardGeneric("predictCoding")
)

setGeneric("refLocsToLocalLocs",
    signature = c("ranges", "txdb", "cdsbytx"),
    function(ranges, txdb, cdsbytx, ...)
    standardGeneric("refLocsToLocalLocs")
)

setGeneric("getTranscriptSeqs", signature = c("query", "subject"),
    function(query, subject, ...)
    standardGeneric("getTranscriptSeqs")
)

### -------------------------------------------------------------------------
### locateVariants 
###

setGeneric("locateVariants", signature = c("query", "subject", "region"),
    function(query, subject, region, ...)
        standardGeneric("locateVariants")
)

### -------------------------------------------------------------------------
### read/write Vcf 
###

setGeneric("readVcf", signature = c("file", "genome", "param"),
    function(file, genome, param, ...)
    standardGeneric("readVcf")
)

setGeneric("readVcfLongForm", signature = c("file", "genome", "param"),
    function(file, genome, param, ...)
    standardGeneric("readVcfLongForm")
)

setGeneric("writeVcf", signature = c("obj", "filename"),
    function(obj, filename, ...)
    standardGeneric("writeVcf")
)

setGeneric("MatrixToSnpMatrix", signature = c("callMatrix", "ref", "alt"),
    function(callMatrix, ref, alt, ...)
    standardGeneric("MatrixToSnpMatrix")
)

### -------------------------------------------------------------------------
### scanVcf 
###

setGeneric("ScanVcfParam", signature = "which",
           function(fixed=character(), info=character(), geno=character(), 
                    trimEmpty=TRUE, which, ...)
           standardGeneric("ScanVcfParam")
)

setGeneric("scanVcfHeader", signature = "file",
           function(file, ...) standardGeneric("scanVcfHeader")
)

setGeneric("scanVcf", signature = c("file", "param"),
           function(file, ..., param) standardGeneric("scanVcf")
)

### -------------------------------------------------------------------------
### VCF class 
###

setGeneric("fixed", signature = "x", 
    function(x) standardGeneric("fixed")
)

setGeneric("fixed<-", signature = c("x", "value"),
    function(x, value) standardGeneric("fixed<-")
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

### -------------------------------------------------------------------------
### VCFHeader class 
###

setGeneric("reference", signature = "x",
    function(x) standardGeneric("reference"),
)

setGeneric("samples", signature = "x",
    function(x) standardGeneric("samples"),
)

setGeneric("header", signature = "x",
    function(x) standardGeneric("header"),
)


setGeneric("contig", signature = "x",
    function(x) standardGeneric("contig"),
)

setGeneric("meta", signature = "x",
    function(x) standardGeneric("meta"),
)

### -------------------------------------------------------------------------
### VAFilter and VAFilterResult 
###

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

