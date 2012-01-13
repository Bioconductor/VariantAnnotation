## VCF 
setClass("VCF",
    contains="SummarizedExperiment",
    representation(
        info="SimpleList"
    ),
)

.valid.VCF <- function(x)
{
    ## FIXME : constraint for Vector, CompressedList, matrix, array
    msg <- NULL
    msg1 <- c("length of info(<VCF>)[[%d]] does not match dim(<VCF>)[1]")

    xlen <- dim(x)[1]
    for (i in seq_len(length(info(x)))) {
        if (class(info(x, i)) == "array")
            len <- dim(info(x, i))[1]
        else
            len <- length(info(x, i))

        if (len != xlen) {
            msg <- c(msg, sprintf(msg1, i))
            next 
        }
    }
    msg 
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


