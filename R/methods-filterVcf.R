setMethod("filterVcf", "character",
    function(file, genome, destination, ..., verbose=TRUE,
             index=FALSE, prefilters=FilterRules(),
             filters=FilterRules(), param=ScanVcfParam())
{
   if (file.exists(destination))
        stop(sprintf("file '%s' exists and will not be over-written", 
             destination))

    tbx <- open(TabixFile(file, yieldSize=100000))
    on.exit(close(tbx))
 

    filterVcf(tbx, genome=genome, destination=destination, ...,
              verbose=verbose, index=index, prefilters=prefilters,
              filters=filters, param=param)
})

.unlistScan <- function(...)
    unlist(scanTabix(...), use.names=FALSE)

.prefilter <-
    function(tbxFile, verbose, prefilters, param, ...)
{
    if (verbose)
        message("starting prefilter")
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
    nTotal <- 0L
    while (length(tbxChunk <- .unlistScan(tbxFile, ..., param=param))) {
        if (verbose)
            message("prefiltering ", nTotal <- nTotal + length(tbxChunk),
                    " records")
        tbxChunk <- subsetByFilter(tbxChunk, prefilters)
        writeLines(tbxChunk, prefiltered)
    }
    close(prefiltered)
    needsClosing <- FALSE

    prefilteredFilename
}

.filter <-
    function(tbxFile, genome, destination, verbose, filters, param, ...)
{
    if (verbose)
        message("starting filter")
    if (!isOpen(tbxFile)) {
        open(tbxFile)
        on.exit(close(tbxFile), add=TRUE)
    }

    filtered <- file(destination, open="w")
    needsClosing <- TRUE
    on.exit(if (needsClosing) close(filtered), add=TRUE)

    nTotal <- 0L
    while (nrow(vcfChunk <- readVcf(tbxFile, genome, ..., param=param))) {
        if (verbose)
            message("filtering ", nTotal <- nTotal + nrow(vcfChunk),
                    " records")
        vcfChunk <- subsetByFilter(vcfChunk, filters)
        writeVcf(vcfChunk, filtered)
    }
    close(filtered)
    needsClosing <- FALSE

    if (verbose)
        message("completed filtering")
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

    if (length(prefilters)) {
        yieldSize <- yieldSize(file)
        file <- .prefilter(file, verbose, prefilters, param, ...)
        if(verbose){
            message(sprintf("prefiltered to %s", file))
            }
        if (length(filters)) {
            ## TabixFile needs to be bgzipped and indexed
            ## FIXME: all records are read at next stage, so no need to index?
            if (verbose)
                message("prefilter compressing and indexing ", sQuote(file))
            gzFilename <- bgzip(file, overwrite = TRUE)
            indexTabix(gzFilename, format = "vcf4")
            file <- TabixFile(gzFilename, yieldSize=yieldSize)
        } else {
            ## file.rename does not work across file systems, so be expensive
            file.copy(file, destination)
            unlink(file)
        }
      }

    if (length(filters)) {
        file <- .filter(file, genome, destination, verbose, filters,
                        param, ...)
        } # if filters
 
    if (index) {
        if (verbose)
            message("compressing and indexing ", sQuote(file))
        gzFilename <- sprintf("%s.gz", destination)
        gzFilename <- bgzip(file, gzFilename, overwrite = TRUE)
        destination <- indexTabix(gzFilename, format = "vcf")
        unlink(file)
    }

    invisible(destination)
})
