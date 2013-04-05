### =========================================================================
### readVcfLongForm methods 
### =========================================================================

## TabixFile

msg <- paste0("'readVcfLongForm' is deprecated. Use 'expand' instead. ",
              "See ?'expand,CollapsedVCF-method'")
setMethod(readVcfLongForm, c(file="TabixFile", genome="character", 
          param="ScanVcfParam"), 
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
    .readVcfLongForm(file, genome, param)
})

setMethod(readVcfLongForm, c(file="TabixFile", genome="character",
          param="GRanges"),
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
    .readVcfLongForm(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcfLongForm, c(file="TabixFile", genome="character",
          param="RangedData"),
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
    .readVcfLongForm(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcfLongForm, c(file="TabixFile", genome="character",
          param="RangesList"),
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
    .readVcfLongForm(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcfLongForm, c(file="TabixFile", genome="character",
          param="missing"), 
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
    .readVcfLongForm(file, genome, param=ScanVcfParam())
})

## character

setMethod(readVcfLongForm, c(file="character", genome="character",
          param="ScanVcfParam"),
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
    file <- .checkTabix(file)
    .readVcfLongForm(file, genome, param)
})

setMethod(readVcfLongForm, c(file="character", genome="character",
          param="missing"),
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
    file <- .checkTabix(file)
    .readVcfLongForm(file, genome, param=ScanVcfParam())
})

setMethod(readVcfLongForm, c(file="character", genome="missing",
          param="missing"),
    function(file, genome, param, ...)
{
    .Deprecated(msg=msg)
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

.readVcfLongForm <- function(file, genome, param = ScanVcfParam(), ...)
{
    .scanVcfToLongForm(scanVcf(file, param=param), file, genome, param)
}

.scanVcfToLongForm <- function(vcf, file, genome, param, ...)
{
    hdr <- scanVcfHeader(file)
    vcf <- .collapseLists(vcf, param)

    ## rowData
    rowData <- vcf[["rowData"]]
    genome(seqinfo(rowData)) <- genome
    values(rowData) <- DataFrame(paramRangeID=vcf[["paramRangeID"]])

    ## fixed fields
    ALT <- .formatALT(vcf[["ALT"]])
    fx <- list(ID=names(rowData), REF=vcf[["REF"]], ALT=ALT, 
               QUAL=vcf[["QUAL"]], FILTER=vcf[["FILTER"]])
    fixed <- DataFrame(fx[lapply(fx, is.null) == FALSE])

    ## info 
    info <- .formatInfo(vcf[["INFO"]], info(hdr))

    ## expand to match ALT
    names(rowData) <-  NULL
    values(rowData) <- append(values(rowData), c(fixed, info))
    rowData <- rep(rowData, elementLengths(ALT))
    values(rowData)[["ALT"]] <- unlist(ALT, use.names=FALSE)
    rowData 
}

.formatGeno <- function(x)
{
    if (length(x) == 0L)
        return(list())
    cls <- lapply(x, class)
    nvar <- dim(x[[1]])[1]
    nsmp <- dim(x[[1]])[2]
    nms <- colnames(x[[1]]) 

    ## collapse matrices and arrays
    for (i in which(cls == "array")) {
        dim(x[[i]]) <- c(nvar*nsmp, dim(x[[i]])[3])
        x[[i]] <- split(x[[i]], seq_len(nvar*nsmp))
    }
    for (i in which(cls == "matrix")) {
        dim(x[[i]]) <- c(nvar*nsmp, 1)
    }
    ## sample names become a data column
    c(list(SAMPLES=rep(nms, each=nvar)), x)
}


