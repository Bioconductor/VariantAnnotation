## VCF 
setClass("VCF",
    contains="SummarizedExperiment",
    representation(
        fixedFields="DataFrame",
        info="DataFrame"
    ),
)

.valid.VCF <- function(x)
{
    xlen <- dim(x)[1]
    if (nrow(values(fixedFields(x))) != xlen)
        stop("'fixedFields' must have the same number of rows as 'rowData'")
    if (nrow(values(info(x))) != xlen)
        stop("'info' must have the same number of rows as 'rowData'")
}

setValidity2("VCF", .valid.VCF)

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


