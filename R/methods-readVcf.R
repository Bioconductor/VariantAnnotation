### =========================================================================
### readVcf methods 
### =========================================================================

setMethod(readVcf, c(file="TabixFile", genome="character", param="ScanVcfParam"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="missing"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="TabixFile", genome="missing", param="missing"), 
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing")
})

setMethod(readVcf, c(file="character", genome="character", param="ScanVcfParam"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="character", genome="character", param="missing"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="character", genome="missing", param="missing"),
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing")
})

.readVcf <- function(file, genome, param = NULL, ...)
{
    if (is.null(param)) 
        vcf <- unpackVcf(scanVcf(file), file)
    else
        vcf <- unpackVcf(scanVcf(file, param=param), file)
    .VcfToSummarizedExperiment(vcf, file, genome)
}

