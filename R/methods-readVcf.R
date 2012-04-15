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
    vcf <- .collapseLists(vcf, param)

    ## rowData
    rowData <- vcf[["rowData"]]
    if (length(rowData))
        genome(seqinfo(rowData)) <- genome 
    values(rowData) <- DataFrame(paramRangeID=vcf[["paramRangeID"]])

    ## fixed fields
    ALT <- .formatALT(vcf[["ALT"]])
    fx <- list(REF=vcf[["REF"]], ALT=ALT, QUAL=vcf[["QUAL"]],
               FILTER=vcf[["FILTER"]])
    fixed <- DataFrame(fx[lapply(fx, is.null) == FALSE]) 

    ## info 
    info <- .formatInfo(vcf[["INFO"]], info(hdr))

    ## colData
    if (length(vcf[["GENO"]]) > 0) {
        samples <- samples(hdr) 
        colData <- DataFrame(Samples=seq_len(length(samples)),
                             row.names=samples)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    VCF(rowData=rowData, colData=colData, exptData=SimpleList(header=hdr),
        fixed=fixed, info=info, geno=SimpleList(vcf[["GENO"]]))
}

