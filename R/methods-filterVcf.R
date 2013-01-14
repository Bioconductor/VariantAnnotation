setMethod("filterVcf", "character",
    function(file, genome, destination, ..., verbose = TRUE,
             index = FALSE, prefilters = FilterRules(),
             filters = FilterRules(), param = ScanVcfParam())
{
    tbx <- open(TabixFile(file, yieldSize=100000))
    on.exit(close(tbx))

    filterVcf(tbx, genome = genome, destination=destination, ...,
              verbose = verbose, index = index, prefilters = prefilters,
              filters = filters, param = param)
})

.unlistScan <- function(...)
    unlist(scanTabix(...), use.names=FALSE)

.prefilter <-
    function(tbxFile, verbose, prefilters, param, ...)
{
    if (!isOpen(tbxFile)) {
        open(tbxFile)
        on.exit(close(tbxFile), add=TRUE)
    }

    prefilteredFilename <- tempfile()
    prefiltered <- file(prefilteredFilename, "w")
    needsClosing <- TRUE
    on.exit(if (needsClosing) close(prefiltered), add=TRUE)
    
    ## copy header
    writeLines(headerTabix(tbxFile)$header, prefiltered)

    ## prefilter
    param <- vcfWhich(param)
    while (length(tbxChunk <- .unlistScan(tbxFile, ..., param=param))) {
        tbxChunk <- subsetByFilter(tbxChunk, prefilters)
        writeLines(tbxChunk, prefiltered)
    }
    close(prefiltered)
    needsClosing <- FALSE

    ## TabixFile needs to be bgzipped and indexed
    ## FIXME: all records are read at next stage, so no need to index?
    gzFilename<- bgzip(prefilteredFilename, overwrite = TRUE)
    indexTabix(gzFilename, format = "vcf")
    unlink(prefilteredFilename)

    TabixFile(gzFilename, yieldSize=yieldSize(tbxFile))
}

.filter <-
    function(tbxFile, genome, destination, verbose, filters, param, ...)
{
    if (!isOpen(tbxFile)) {
        open(tbxFile)
        on.exit(close(tbxFile), add=TRUE)
    }

    filtered <- file(destination, open="a")
    needsClosing <- TRUE
    on.exit(if (needsClosing) close(filtered), add=TRUE)

    while (nrow(vcfChunk <- readVcf(tbxFile, genome, ..., param=param))) {
        vcfChunk <- subsetByFilter(vcfChunk, filters)
        writeVcf(vcfChunk, filtered)
    }
    close(filtered)
    needsClosing <- FALSE

    destination
}

setMethod("filterVcf", "TabixFile",
    function(file, genome, destination, ..., verbose = TRUE,
             index = FALSE,
             prefilters = FilterRules(), filters = FilterRules(),
             param = ScanVcfParam())
{
    if (!isSingleString(destination))
        stop("'destination' must be character(1)")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (!isTRUEorFALSE(index))
        stop("'index' must be TRUE or FALSE")

    if (!length(prefilters) && !length(filters))
        stop("no 'prefilters' or 'filters' specified")

    if (length(prefilters))
        file <- .prefilter(file, verbose, prefilters, param, ...)

    if (length(filters))
        file <- .filter(file, genome, destination, verbose, filters,
                        param, ...)

    if (index) {
        gzFilename <- sprintf("%s.gz", destination)
        gzFilename <- bgzip(file, gzFilename, overwrite = TRUE)
        destination <- indexTabix(gzFilename, format = "vcf")
        unlink(file)
    } else {
        if (file != destination)
            file.rename(file, destination)
    }

    invisible(destination)
})

## Deprecated

dbSNPFilter <-
    function(dbSNP=character(0), .name="dbSNPFilter")
{
    .Deprecated("filterVcf")
}

regionFilter <-
    function(txdb, region="coding", .name="regionFilter")
{
    .Deprecated("filterVcf")
}
