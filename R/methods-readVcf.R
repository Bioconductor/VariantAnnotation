### =========================================================================
### readVcf methods 
### =========================================================================

## TabixFile

setMethod(readVcf, c(file="TabixFile", genome="character", 
          param="ScanVcfParam"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="GRanges"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="RangedData"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="GRangesList"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="RangesList"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="missing"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam())
})

## character

setMethod(readVcf, c(file="character", genome="character",
          param="ScanVcfParam"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="character", genome="character",
          param="missing"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome, param=ScanVcfParam())
})

setMethod(readVcf, c(file="character", genome="missing",
          param="missing"),
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing") 
})

.checkTabix <- function(x)
{
    if (1L != length(x)) 
        stop("'x' must be character(1)")
    if (grepl("\\.tbi$", x))
        TabixFile(sub("\\.tbi", "", x))
    else 
        x 
}

.readVcf <- function(file, genome, param, ...)
{
    .scanVcfToVCF(scanVcf(file, param=param), file, genome, param)
}

.scanVcfToVCF <- function(vcf, file, genome, param, ...)
{
    hdr <- scanVcfHeader(file)
    if (length(vcf[[1]]$GENO) > 0L)
        colnms <- colnames(vcf[[1]]$GENO[[1]])
    else
        colnms <- NULL
    vcf <- .collapseLists(vcf, param)

    ## rowData
    rowData <- vcf$rowData
    if (length(rowData))
        genome(seqinfo(rowData)) <- genome 
    values(rowData) <- DataFrame(vcf["paramRangeID"])

    ## fixed fields
    fx <- vcf[c("REF", "ALT", "QUAL", "FILTER")]
    fx$ALT <- .formatALT(fx$ALT)
    fixed <- DataFrame(fx[!sapply(fx, is.null)]) 

    ## info 
    info <- .formatInfo(vcf$INFO, info(hdr), length(rowData))

    ## colData
    colData <- DataFrame(Samples=seq_along(colnms), row.names=colnms)

    VCF(rowData=rowData, colData=colData, exptData=SimpleList(header=hdr),
        fixed=fixed, info=info, geno=SimpleList(vcf$GENO))
}

