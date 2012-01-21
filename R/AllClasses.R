## VCF 

setClass("VCF",
    contains="SummarizedExperiment",
    representation(
        fixedFields="DataFrame",
        info="DataFrame"
    ),
)

.valid.VCF.fixedFields <- function(x)
{
    xlen <- dim(x)[1]
    ff <- values(fixedFields(x))
    if (nrow(ff) != xlen)
        return(paste("'fixedFields(x)' and 'rowData(x) must have the same ",
               "number of rows", sep=""))

    if (!identical(DataFrame(), ff)) {
        if (!all(names(ff) %in% c("REF", "ALT", "QUAL", "FILTER")))
            return(paste("'values(fixedFields(x))' colnames must be 'REF', ",
                   "'ALT', 'QUAL' and 'FILTER'", sep=""))
        if (!is(ff$REF, "DNAStringSet"))
            return("'values(fixedFields(x))[['REF']] must be a DNAStringSet")
        if (!is(ff$ALT, "DNAStringSetList") && !is(ff$ALT, "CharacterList"))
            return(paste("'values(fixedFields(x))[['ALT']] must be a ",
                   "DNAStringSetList or a CharacterList", sep=""))
        if (!is(ff$QUAL, "numeric"))
            return("'values(fixedFields(x))[['QUAL']] must be numeric")
        if (!is(ff$FILTER, "character"))
            return("'values(fixedFields(x))[['FILTER']] must be a character")
    }
    NULL
}

.valid.VCF.info <- function(x)
{
    xlen <- dim(x)[1]
    if (nrow(values(info(x))) != xlen)
        return("'info' must have the same number of rows as 'rowData'")
    NULL
}

.valid.VCF <- function(x)
{
    c(.valid.VCF.fixedFields(x),
      .valid.VCF.info(x))
}

setValidity2("VCF", .valid.VCF, where=asNamespace("VariantAnnotation"))

setClass("ScanVcfParam", contains="ScanBVcfParam")


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


