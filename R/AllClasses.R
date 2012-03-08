## VCF 

setClass("VCF",
    contains="SummarizedExperiment",
    representation(
        fixed="DataFrame",
        info="DataFrame"
    ),
)

.valid.VCF.fixed <- function(x)
{
    xlen <- dim(x)[1]
    ffld <- values(fixed(x))
    nms <- names(ffld)

    if (length(ffld) != 0) {
        if (nrow(ffld) != xlen)
            return(paste("'fixed(x)' and 'rowData(x) must have the same ",
                   "number of rows", sep=""))
        if (!all(nms %in% c("paramRangeID", "REF", "ALT", "QUAL", "FILTER")))
            return(paste("'values(fixed(x))' colnames must be ",
                   "'REF', 'ALT', 'QUAL' and 'FILTER'", sep=""))
        if ("REF" %in% nms) 
            if (!is(ffld$REF, "DNAStringSet"))
                return("'values(fixed(x))[['REF']] must be a DNAStringSet")
        if ("ALT" %in% nms) 
            if (!is(ffld$ALT, "DNAStringSetList") && !is(ffld$ALT, "CharacterList"))
                return(paste("'values(fixed(x))[['ALT']] must be a ",
                       "DNAStringSetList or a CharacterList", sep=""))
        if ("QUAL" %in% nms) 
            if (!is(ffld$QUAL, "numeric"))
                return("'values(fixed(x))[['QUAL']] must be numeric")
        if ("FILTER" %in% nms) 
            if (!is(ffld$FILTER, "character"))
                return("'values(fixed(x))[['FILTER']] must be a character")
    }
    NULL
}

.valid.VCF.info <- function(x)
{
    xlen <- dim(x)[1]
    info <- values(info(x))
    if (length(info) != 0)
        if (nrow(values(info(x))) != xlen)
            return("'info' must have the same number of rows as 'rowData'")
    NULL
}

.valid.VCF <- function(x)
{
    c(.valid.VCF.fixed(x),
      .valid.VCF.info(x))
}

setValidity2("VCF", .valid.VCF, where=asNamespace("VariantAnnotation"))

setClass("ScanVcfParam", contains="ScanBVcfParam")

## vcfHeader 

setClass("VCFHeader", 
    representation(
        reference="character",
        samples="character",
        header="SimpleDataFrameList"
    )
)

## VAFilters 

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

## SIFT and PolyPhen 

.SIFTDb <- setRefClass("SIFTDb",
    contains = "AnnotationDb",
)

.PolyPhenDb <- setRefClass("PolyPhenDb",
    contains = "AnnotationDb",
)


