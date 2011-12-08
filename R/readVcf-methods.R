### =========================================================================
### readVcf methods 
### =========================================================================

setMethod(readVcf, c(file="TabixFile", genome="character", param="RangesList"),
    function(file, genome, ..., param)
{
    .readVcf(file, genome, param=param)
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="RangedData"),
    function(file, genome, ..., param)
{
    .readVcf(file, genome, param=param)
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="GRanges"),
    function(file, genome, ..., param)
{
    .readVcf(file, genome, param=param)
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="ScanVcfParam"), 
    function(file, genome, ..., param)
{
    if (length(vcfWhich(param)) == 0)
        file <- path(file)
    .readVcf(file, genome, param=param)
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="missing"), 
    function(file, genome, ..., param)
{
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="character", genome="character", param="missing"),
    function(file, genome, ..., param)
{
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="character", genome="character", param="ANY"),
    function(file, genome, ..., param)
{
    .readVcf(file, genome, param=param)
})

.readVcf <- function(file, genome, ..., param = NULL)
{
    if (is.null(param)) 
        vcf <- unpackVcf(scanVcf(file), file)
    else
        vcf <- unpackVcf(scanVcf(file, param=param), file)
    .VcfToSummarizedExperiment(vcf, file, genome, param=param)
}

