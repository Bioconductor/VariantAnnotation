setMethod("filterVcf", "character",
    function(file, genome, destination, ..., verbose = TRUE,
             index = FALSE, filters = FilterRules(),
             param = ScanVcfParam())
{
    tbx <- open(TabixFile(file, yieldSize=10000))
    on.exit(close(tbx))

    filterVcf(tbx, genome = genome, destination=destination, ...,
              verbose = verbose, index = index, filters = filters,
              param = param)
})

setMethod("filterVcf", "TabixFile",
    function(file, genome, destination, ..., verbose = TRUE,
             index = FALSE, filters = FilterRules(),
             param = ScanVcfParam())
{
    if (!isSingleString(destination))
        stop("'destination' must be character(1)")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (!isTRUEorFALSE(index))
        stop("'index' must be TRUE or FALSE")

    ## desination: character(1) file path
    con <- file(destination, open="a")
    needsClosing <- TRUE
    on.exit(if (needsClosing) close(con))

    if (!isOpen(file)) {
        open(file)
        on.exit(close(file), add=TRUE)
    }
    if (verbose) {
        wd <- options()$width
        pbi <- 0L
        pb <- txtProgressBar(max=wd)
    }
    while (nrow(vcfChunk <- readVcf(file, genome, ..., param=param))) {
        if (verbose) {
            pbi <- (pbi + 1L) %% wd
            if (pbi == 0)
                message()
            setTxtProgressBar(pb, pbi)
        }
        vcfChunk <- subsetByFilter(vcfChunk, filters)
        writeVcf(vcfChunk, con)
    }
    if (verbose)
        message()
    close(con)                          # so bgzip is on a closed file
    needsClosing <- FALSE

    if (index) {
        filenameGZ <- bgzip(destination, overwrite = TRUE)
        indexTabix(filenameGZ, format = "vcf")
        unlink(destination)
        invisible(filenameGZ)
    } else {
        invisible(destination)
    }
})
