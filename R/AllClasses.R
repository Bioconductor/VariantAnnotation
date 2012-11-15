### ------------------------------------------------------------------------- 
### VCF class hierarchy 
###

setClass("VCF",
    contains=c("VIRTUAL", "SummarizedExperiment"),
    representation(
        fixed="DataFrame",
        info="DataFrame"
    )
)

setClass("CollapsedVCF", contains="VCF") ## ALT is DNAStrinsSetList

setClass("ExpandedVCF", contains="VCF")  ## ALT is DNAStringSet 

### Coercion: 
### Recursion problem in an automatically generated coerce method requires
### that we handle coercion from subclasses to SummarizedExperiment.

setAs("ExpandedVCF", "SummarizedExperiment", 
    def = function(from)
    { 
        if (strict) {
            force(from)
            class(from) <- "SummarizedExperiment"
        }
        from
    },
    replace = function(from, to, value)
    {
        firstTime <- TRUE
        for (nm in slotNames(value)) {
            v <- slot(value, nm)
            if (firstTime) {
                slot(from, nm, FALSE) <- v
                firstTime <- FALSE
            } else {
                `slot<-`(from, nm, FALSE, v)
            }
        }
        from
    }
)

setAs("CollapsedVCF", "SummarizedExperiment", 
    def = function(from)
    { 
        if (strict) {
            force(from)
            class(from) <- "SummarizedExperiment"
        }
        from
    },
    replace = function(from, to, value)
    {
        firstTime <- TRUE
        for (nm in slotNames(value)) {
            v <- slot(value, nm)
            if (firstTime) {
                slot(from, nm, FALSE) <- v
                firstTime <- FALSE
            } else {
                `slot<-`(from, nm, FALSE, v)
            }
        }
        from
    }
)

### Validity 

.valid.VCF.fixed <- function(object)
{
    xlen <- dim(object)[1]
    ffld <- slot(object, "fixed")
    nms <- names(ffld)

    if (nrow(ffld) != 0) {
        if (nrow(ffld) != xlen)
            return(paste("'fixed(object)' and 'rowData(object) must have the same ",
                   "number of rows", sep=""))
        if (!all(nms %in% c("paramRangeID", "REF", "ALT", "QUAL", "FILTER")))
            return(paste("'mcols(fixed(object))' colnames must be ",
                   "'REF', 'ALT', 'QUAL' and 'FILTER'", sep=""))
        if ("REF" %in% nms) 
            if (!is(ffld$REF, "DNAStringSet"))
                return("'ref(object)' must be a DNAStringSet")
        if ("QUAL" %in% nms) 
            if (!is(ffld$QUAL, "numeric"))
                return("'qual(object)' must be numeric")
        if ("FILTER" %in% nms) 
            if (!is(ffld$FILTER, "character"))
                return("'filt(object)' must be a character")
    }
    NULL
}

.valid.VCF.info <- function(object)
{
    xlen <- nrow(object)
    i <- slot(object, "info")
    if (length(i) != 0)
        if (nrow(i) != xlen)
            return("'info' must have the same number of rows as 'rowData'")
    NULL
}

.valid.CollapsedVCF.alt <- function(object)
{
    ffld <- slot(object, "fixed")
    if (length(ffld) != 0L) {
        alt <- alt(object)
        if (length(alt) != 0L)
            if (!is(alt, "DNAStringSetList") && !is(alt, "CharacterList")) 
                return(paste("'alt(object)' must be a DNAStringSetList or a ",
                       "CharacterList", sep=""))
    } 
    NULL
}
 
.valid.ExpandedVCF.alt <- function(object)
{
    ffld <- slot(object, "fixed")
    if (length(ffld) != 0L) {
        alt <- alt(object)
        if (length(alt) != 0L)
            if (!is(alt, "DNAStringSet") && !is(alt, "CharacterList")) 
                return(paste("'alt(object)' must be a DNAStringSet or a ",
                       "CharacterList", sep=""))
    } 
    NULL
}

.valid.CollapsedVCF <- function(object)
{
    c(.valid.VCF.fixed(object),
      .valid.VCF.info(object),
      .valid.CollapsedVCF.alt(object))
}

.valid.ExpandedVCF <- function(object)
{
    c(.valid.VCF.fixed(object),
      .valid.VCF.info(object),
      .valid.ExpandedVCF.alt(object))
}
 
setValidity("CollapsedVCF", .valid.CollapsedVCF, where=asNamespace("VariantAnnotation"))
setValidity("ExpandedVCF", .valid.ExpandedVCF, where=asNamespace("VariantAnnotation"))

### ------------------------------------------------------------------------- 
### ScanVcfParam class 
###

setClass("ScanVcfParam", contains="ScanBVcfParam")

### ------------------------------------------------------------------------- 
### VCFHeader class 
###

setClass("VCFHeader", 
    representation(
        reference="character",
        samples="character",
        header="SimpleDataFrameList"
    )
)

### ------------------------------------------------------------------------- 
### SIFT and PolyPhen classes 
###

.SIFTDb <- setRefClass("SIFTDb", contains = "AnnotationDb")

.PolyPhenDb <- setRefClass("PolyPhenDb", contains = "AnnotationDb")

### ------------------------------------------------------------------------- 
### VariantType class hierarchy 
###

setClass("VariantType",
    representation(
        "VIRTUAL"
    )
)

setMethod("show", "VariantType",
    function(object) 
    {
        cat("class:", class(object), "\n")
    }
)

setClass("CodingVariants", contains="VariantType")

setClass("IntronVariants", contains="VariantType")

setClass("IntergenicVariants", contains="VariantType")

setClass("ThreeUTRVariants", contains="VariantType")

setClass("FiveUTRVariants", contains="VariantType")

setClass("SpliceSiteVariants", contains="VariantType")

setClass("PromoterVariants", 
    contains="VariantType",
    representation(upstream="integer",
                   downstream="integer")
)

.promoterValidity <- function(object, ...)
{
    if (any((upstream(object) < 0)  | (downstream(object) < 0)))
        return("'upstream' and 'downstream' must be integers >= 0")
    NULL 
}

setValidity("PromoterVariants",
    function(object)
        .promoterValidity(object)
)

setClass("AllVariants", 
    contains="VariantType",
    representation(upstream="integer",
                   downstream="integer")
)

setValidity("AllVariants",
    function(object)
        .promoterValidity(object)
)
