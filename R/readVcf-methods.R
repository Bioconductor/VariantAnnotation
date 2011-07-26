## TODO :
## imprecise variants ALT field may contain descriptive alleles
## imprecise should be marked by imprecise flag in INFO field
## strand bias in INFO$sb

## FIXME : if user does not supply param, either 
##         (1) make a default param for scanTabix such that
##             all records are brought over
##         (2) use scanBcf(file, character(0))

setMethod("readVcf", c(file="character", param="missing"),
    function(file, ..., param, raw=FALSE)
{
    vcf <- scanBcf(file, character(0))
    .VcfToSummarizedExperiment(vcf, raw=raw)
})

setMethod("readVcf", c(file="character", param="ANY"),
    function(file, index=paste(file, "tbi", sep="."), ..., param, raw=FALSE)
{
    tf <- TabixFile(file, index=paste(file, "tbi", sep="."))
    callGeneric(tf, param=param, raw=raw)
})

setMethod("readVcf", c(file="TabixFile", param="GRanges"),
    function(file, ..., param, raw=FALSE)
{
    .readVcf(file, param=param, raw=raw)
})

setMethod("readVcf", c(file="TabixFile", param="RangedData"),
    function(file, ..., param, raw=FALSE)
{
    .readVcf(file, param=param, raw=raw)
})

setMethod("readVcf", c(file="TabixFile", param="RangesList"),
    function(file, ..., param, raw=FALSE)
{
    .readVcf(file, param=param, raw=raw)
})


.readVcf <- function(file, param, raw, ...)
{
    tbx <- scanTabix(file, param=param)
    vcf <- .parseTabix(tbx, param=param)
    .VcfToSummarizedExperiment(vcf, raw=raw)
}
