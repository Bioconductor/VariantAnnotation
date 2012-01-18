### =========================================================================
### readVcf methods 
### =========================================================================

## TabixFile

setMethod(readVcf, c(file="TabixFile", genome="character", param="ScanVcfParam"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="GRanges"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="RangedData"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="RangesList"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="missing"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="TabixFile", genome="missing", param="ANY"), 
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing")
})


## character

setMethod(readVcf, c(file="character", genome="character", param="ScanVcfParam"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="character", genome="character", param="GRanges"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="character", genome="character", param="RangedData"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="character", genome="character", param="RangesList"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="character", genome="character", param="missing"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="character", genome="missing", param="ANY"),
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing")
})


## .readVcf internal

.readVcf <- function(file, genome, param, ...)
{
    if (missing(param)) {
        .scanVcfToVCF(scanVcf(file), file, genome)
    } else {
        if (!identical(character(), vcfAsGRanges(param)))
            .scanVcfToLongGRanges(scanVcf(file, param=param), file, genome, param)
        else
            .scanVcfToVCF(scanVcf(file, param=param), file, genome)
    }
}

