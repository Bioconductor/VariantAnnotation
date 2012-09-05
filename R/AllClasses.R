### ------------------------------------------------------------------------- 
### VCF class 
###

setClass("VCF",
    contains="SummarizedExperiment",
    representation(
        fixed="DataFrame",
        info="DataFrame"
    )
)

.valid.VCF.fixed <- function(object)
{
    xlen <- dim(object)[1]
    ffld <- slot(object, "fixed")
    nms <- names(ffld)

    if (length(ffld) != 0) {
        if (nrow(ffld) != xlen)
            return(paste("'fixed(object)' and 'rowData(object) must have the same ",
                   "number of rows", sep=""))
        if (!all(nms %in% c("paramRangeID", "REF", "ALT", "QUAL", "FILTER")))
            return(paste("'values(fixed(object))' colnames must be ",
                   "'REF', 'ALT', 'QUAL' and 'FILTER'", sep=""))
        if ("REF" %in% nms) 
            if (!is(ffld$REF, "DNAStringSet"))
                return("'values(fixed(object))[['REF']] must be a DNAStringSet")
        if ("ALT" %in% nms) 
            if (!is(ffld$ALT, "DNAStringSetList") && !is(ffld$ALT, "CharacterList"))
                return(paste("'values(fixed(object))[['ALT']] must be a ",
                       "DNAStringSetList or a CharacterList", sep=""))
        if ("QUAL" %in% nms) 
            if (!is(ffld$QUAL, "numeric"))
                return("'values(fixed(object))[['QUAL']] must be numeric")
        if ("FILTER" %in% nms) 
            if (!is(ffld$FILTER, "character"))
                return("'values(fixed(object))[['FILTER']] must be a character")
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

.valid.VCF <- function(object)
{
    c(.valid.VCF.fixed(object),
      .valid.VCF.info(object))
}

setValidity("VCF", .valid.VCF, where=asNamespace("VariantAnnotation"))

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
### VAFilters classes 
###

setClass(".VAUtil",
    representation("VIRTUAL")
)

.vaValidity <- function(object) TRUE
setClass("VAFilter",
    contains=c("function", ".VAUtil"),
    representation(
        name="ScalarCharacter"
    ),
    validity=.vaValidity
)

setClass("VAFilterResult",
    contains=c("logical", ".VAUtil"),
    representation(
      name="ScalarCharacter",
      stats="data.frame"
    )
)

### ------------------------------------------------------------------------- 
### SIFT and PolyPhen classes 
###

.SIFTDb <- setRefClass("SIFTDb",
    contains = "AnnotationDb",
)

.PolyPhenDb <- setRefClass("PolyPhenDb",
    contains = "AnnotationDb",
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
