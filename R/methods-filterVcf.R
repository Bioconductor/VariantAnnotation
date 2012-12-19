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
    function(tbxFile, verbose, index, prefilters, param, ...)
{
    if (!isOpen(tbxFile)) {
        open(tbxFile)
        on.exit(close(tbxFile), add=TRUE)
    }

    prefilteredFilename <- tempfile()
    prefiltered <- file(prefilteredFilename, "w")
    needsClosing <- TRUE
    on.exit(if (needsClosing) close(prefiltered))
    
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

    if (index) {
        prefilteredFilename <- bgzip(prefilteredFilename, overwrite = TRUE)
        indexTabix(prefilteredFilename, format = "vcf")
    }

    TabixFile(prefilteredFilename, yieldSize=yieldSize(tbxFile))
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
    on.exit(if (needsClosing) close(filtered))

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

    ## prefilters
    if (length(prefilters)) {
        doIndex <- index || (length(filters) != 0)
        file <- .prefilter(file, verbose, doIndex, prefilters, param,
                           ...)
    }
    
    ## filters
    if (length(filters))
        file <- .filter(file, genome, destination, verbose, filters,
                        param, ...)

    ## desination: character(1) file path
    if (index) {
        filenameGZ <- bgzip(file, overwrite = TRUE)
        indexTabix(filenameGZ, format = "vcf")
        unlink(file)
        invisible(filenameGZ)
    } else {
        invisible(file)
    }
})
