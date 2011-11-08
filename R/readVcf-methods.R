### =========================================================================
### readVcf methods 
### =========================================================================

setMethod(readVcf, c(file="TabixFile", param="RangesList"),
    function(file, ..., param)
{
    .readVcf(file, param=param)
})

setMethod(readVcf, c(file="TabixFile", param="RangedData"),
    function(file, ..., param)
{
    .readVcf(file, param=param)
})

setMethod(readVcf, c(file="TabixFile", param="GRanges"),
    function(file, ..., param)
{
    .readVcf(file, param=param)
})

setMethod(readVcf, c(file="TabixFile", param="ScanVcfParam"), 
    function(file, ..., param)
{
    .readVcf(file, param=param)
})

setMethod(readVcf, c(file="TabixFile", param="missing"), 
    function(file, ..., param)
{
    .readVcf(file)
})

setMethod(readVcf, c(file="character", param="missing"),
    function(file, ..., param)
{
    .readVcf(file)
})

setMethod(readVcf, c(file="character", param="ANY"),
    function(file, ..., param)
{
    .readVcf(file, param=param)
})

.readVcf <- function(file, ..., param)
{
    if (missing(param))
        vcf <- unpackVcf(scanVcf(file), file)
    else
        vcf <- unpackVcf(scanVcf(file, param=param), file)
    .VcfToSummarizedExperiment(vcf, file)
}

