### -------------------------------------------------------------------------
### predictCoding 
###

setGeneric("predictCoding", 
    signature=c("query", "subject", "seqSource", "varAllele"),
    function(query, subject, seqSource, varAllele, ...)
        standardGeneric("predictCoding")
)

setGeneric("getTranscriptSeqs", signature=c("query", "subject"),
    function(query, subject, ...)
        standardGeneric("getTranscriptSeqs")
)

### -------------------------------------------------------------------------
### locateVariants 
###

setGeneric("locateVariants", signature=c("query", "subject", "region"),
    function(query, subject, region, ...)
        standardGeneric("locateVariants")
)

### -------------------------------------------------------------------------
### summarizeVariants 
###

setGeneric("summarizeVariants", signature=c("query", "subject", "mode"),
    function(query, subject, mode, ...)
        standardGeneric("summarizeVariants")
)

### -------------------------------------------------------------------------
### VariantRegion classes 
###

setGeneric("upstream", signature="x",
    function(x) standardGeneric("upstream")
)

setGeneric("upstream<-", signature="x",
    function(x, value) standardGeneric("upstream<-")
)

setGeneric("downstream", signature="x",
    function(x) standardGeneric("downstream")
)

setGeneric("downstream<-", signature="x",
    function(x, value) standardGeneric("downstream<-")
)

setGeneric("promoter", signature="x",
    function(x) standardGeneric("promoter")
)

setGeneric("promoter<-", signature="x",
    function(x, value) standardGeneric("promoter<-")
)

setGeneric("intergenic", signature="x",
    function(x) standardGeneric("intergenic")
)

setGeneric("intergenic<-", signature="x",
    function(x, value) standardGeneric("intergenic<-")
)

setGeneric("idType", signature="x",
    function(x) standardGeneric("idType")
)

setGeneric("idType<-", signature="x",
    function(x, value) standardGeneric("idType<-")
)

### -------------------------------------------------------------------------
### read/write Vcf 
###

setGeneric("readVcf", signature=c("file", "param"),
    function(file, genome, param, ...)
    standardGeneric("readVcf")
)

setGeneric("writeVcf", signature=c("obj", "filename"),
    function(obj, filename, ...)
    standardGeneric("writeVcf")
)


### -------------------------------------------------------------------------
### scanVcf 
###

setGeneric("ScanVcfParam", signature="which",
           function(fixed=character(), info=character(), geno=character(), 
                    samples=character(), trimEmpty=TRUE, which, ...)
           standardGeneric("ScanVcfParam")
)

setGeneric("scanVcfHeader", signature="file",
           function(file, ...) standardGeneric("scanVcfHeader")
)

setGeneric("scanVcf", signature=c("file", "param"),
           function(file, ..., param) 
               standardGeneric("scanVcf")
)

### -------------------------------------------------------------------------
### filterVcf
###

setGeneric("filterVcf", signature="file",
           function(file, genome, destination, ..., verbose=FALSE,
                    index=FALSE, prefilters=FilterRules(),
                    filters=FilterRules(), param=ScanVcfParam())
           standardGeneric("filterVcf")
)

### -------------------------------------------------------------------------
### VCF class 
###

setGeneric("fixed", signature="x", 
    function(x) standardGeneric("fixed")
)

setGeneric("fixed<-", signature=c("x", "value"),
    function(x, value) standardGeneric("fixed<-")
)

setGeneric("ref", signature="x", 
    function(x) standardGeneric("ref")
)

setGeneric("ref<-", signature=c("x", "value"),
    function(x, value) standardGeneric("ref<-")
)

setGeneric("alt", signature="x", 
    function(x) standardGeneric("alt")
)

setGeneric("alt<-", signature=c("x", "value"),
    function(x, value) standardGeneric("alt<-")
)

setGeneric("qual", signature="x", 
    function(x) standardGeneric("qual")
)

setGeneric("qual<-", signature=c("x", "value"),
    function(x, value) standardGeneric("qual<-")
)

setGeneric("filt", signature="x", 
    function(x) standardGeneric("filt")
)

setGeneric("filt<-", signature=c("x", "value"),
    function(x, value) standardGeneric("filt<-")
)

setGeneric("info", signature="x", 
    function(x, ..., row.names = TRUE) standardGeneric("info")
)

setGeneric("info<-", signature=c("x", "value"),
    function(x, value) standardGeneric("info<-")
)

setGeneric("geno", signature=c("x", "i"),
    function(x, i, ..., withDimnames=TRUE) 
    standardGeneric("geno"),
)

setGeneric("geno<-", signature=c("x", "i", "value"),
    function(x, i, ..., value) 
    standardGeneric("geno<-")
)

### -------------------------------------------------------------------------
### VCFHeader class 
###

setGeneric("reference", signature="x",
    function(x) standardGeneric("reference"),
)

setGeneric("header", signature="x",
    function(x) standardGeneric("header"),
)

setGeneric("header<-", signature=c("x", "value"),
    function(x, value) standardGeneric("header<-"),
)

setGeneric("contig", signature="x",
    function(x) standardGeneric("contig"),
)

setGeneric("meta", signature="x",
    function(x) standardGeneric("meta"),
)

setGeneric("meta<-", signature=c("x", "value"),
    function(x, value) standardGeneric("meta<-"),
)

### -------------------------------------------------------------------------
### snp encoding methods 
###

setGeneric("genotypeToSnpMatrix", signature="x",
    function(x, ...)
    standardGeneric("genotypeToSnpMatrix")
)

setGeneric("snpSummary", function(x, ...) standardGeneric("snpSummary") )

### -------------------------------------------------------------------------
### isSNV helpers 
###

setGeneric("isSNV", signature="x",
    function(x, ...)
    standardGeneric("isSNV")
)

setGeneric("isInsertion", signature="x",
    function(x, ...)
    standardGeneric("isInsertion")
)

setGeneric("isDeletion", signature="x",
    function(x, ...)
    standardGeneric("isDeletion")
)

setGeneric("isIndel", signature="x",
    function(x, ...)
    standardGeneric("isIndel")
)

setGeneric("isDelins", signature="x",
    function(x, ...)
    standardGeneric("isDelins")
)

setGeneric("isTransition", signature="x",
    function(x, ...)
    standardGeneric("isTransition")
)

setGeneric("isSV", signature="x",
    function(x, ...)
    standardGeneric("isSV")
)

setGeneric("isSVPrecise", signature="x",
    function(x, ...)
    standardGeneric("isSVPrecise")
)

setGeneric("isSubstitution", signature="x",
    function(x, ...)
    standardGeneric("isSubstitution")
)

### -------------------------------------------------------------------------
### VRanges class 
###

setGeneric("totalDepth", function(x, ...) standardGeneric("totalDepth"))
setGeneric("altDepth", function(x, ...) standardGeneric("altDepth"))
setGeneric("refDepth", function(x, ...) standardGeneric("refDepth"))
setGeneric("softFilterMatrix",
           function(x, value) standardGeneric("softFilterMatrix"))
setGeneric("softFilterMatrix<-",
           function(x, value) standardGeneric("softFilterMatrix<-"))
setGeneric("hardFilters",
           function(x, value) standardGeneric("hardFilters"))
setGeneric("hardFilters<-",
           function(x, value) standardGeneric("hardFilters<-"))
setGeneric("called", function(x, ...) standardGeneric("called"))
setGeneric("altFraction", function(x, ...) standardGeneric("altFraction"))
setGeneric("refFraction", function(x, ...) standardGeneric("refFraction"))
setGeneric("asVCF", function(x, ...) standardGeneric("asVCF"))
setGeneric("tabulate", signature = "bin", # BiocGenerics?
           function(bin, nbins = max(1L, bin, na.rm = TRUE))
           standardGeneric("tabulate"))

### -------------------------------------------------------------------------
### VRangesList class 
###

setGeneric("stackSamples", function(x, ...) standardGeneric("stackSamples"))
