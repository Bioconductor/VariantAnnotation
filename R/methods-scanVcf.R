### =========================================================================
### scanVcf methods 
### =========================================================================

.vcf_usertag <-
    function(map, tag, nm, ...)
{
    if (!identical(character(), tag))
        if (1L == length(tag) && is.na(tag)) {
            map[] <- list(list("0", NULL))
        } else {
            map[!names(map) %in% tag] <- list(list("0", NULL))
            if (!all(tag %in% names(map)))
                warning("ScanVcfParam '", nm, "' fields not present: '",
                        paste(tag[!tag %in% names(map)], collapse="' ' "),
                        "'")
        }
    map
}

.vcf_fixed <-
    function(tag, ...)
{
    map <- list(ALT=list("A", character()), QUAL=list("1", numeric()),
                FILTER=list("1", character()))
    c(list(rowData=NULL, REF=NULL), .vcf_usertag(map, tag, "fixed"))
}

.vcf_map <-
    function(fmt, tag, ...)
{
    map <- Map(function(n, t) {
        if (t == "Flag") n <- "1"
        t <- switch(t, String=character(), Integer=integer(),
               Float=numeric(), Flag=logical() )
        list(n, t)
    }, fmt$Number, fmt$Type)
    names(map) <- rownames(fmt)

    ## user selected
    .vcf_usertag(map, tag, ...)
}

.vcf_scan_header_maps <-
    function(file, fixed, info, geno, samples)
{
    hdr <- suppressWarnings(scanVcfHeader(file))
    samp <- samples(hdr)
    if (identical(character(), samples)) {
        smap <- !logical(length(samp))
    } else if (isTRUE(is.na(samples))) {
        smap <- logical(length(samp))
    } else {
        if (length(nx <- samples[!samples %in% samp]) > 0L)
            warning(paste("samples", sQuote(nx), "not found in file"))
        smap <- samp %in% samples
    }
    names(smap) <- samp
    fmap <- .vcf_fixed(fixed)
    imap <- .vcf_map(info(hdr), info, nm="info")
    if (0L == length(imap))
        imap <- list(list("1", character()))
    gmap <- .vcf_map(geno(hdr), geno, nm="geno")
    list(hdr=hdr, samples=smap, fmap=fmap, imap=imap, gmap=gmap)
}

.vcf_scan <- 
   function(file, ..., fixed=character(), info=character(),
            geno=character(), samples=character(), param)
{
    tryCatch({
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }
        maps <- .vcf_scan_header_maps(file, fixed, info, geno, samples)
        tbxstate <- maps[c("samples", "fmap", "imap", "gmap")]
        tbxsym <- getNativeSymbolInfo(".tabix_as_vcf",
                                      "VariantAnnotation")
        scanTabix(file, ..., param=param, tbxsym=tbxsym, tbxstate=tbxstate)
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             path(file), call. = FALSE)
    })
}

.vcf_scan_character <-
    function(file, ..., fixed=character(), info=character(),
             geno=character(), samples=character(), 
             yieldSize=100000L)
{
    tryCatch({
        file <- normalizePath(path.expand(file))
        if (!file.exists(file))
            stop("file does not exist")
        maps <- .vcf_scan_header_maps(file, fixed, info, geno, samples)
        result <- .Call(.scan_vcf_character, file,
                        as.integer(yieldSize), maps$samples,
                        maps$fmap, maps$imap, maps$gmap)
        setNames(result, "*:*-*")
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             file, call. = FALSE)
    })
}

.vcf_scan_connection <-
    function(file, ..., fixed=character(), info=character(),
             geno=character(), samples=character())
{
    tryCatch({
        file <- summary(file)$description
        maps <- .vcf_scan_header_maps(file, fixed, info, geno, samples)

        txt <- readLines(file, ...)
        txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
        result <- .Call(.scan_vcf_connection, txt, maps$samples,
                        maps$fmap, maps$imap, maps$gmap)
        setNames(result, "*:*-*")
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             summary(file)$description, call. = FALSE)
    })
}

setMethod(scanVcf, c("TabixFile", "RangesList"),
    function(file, ..., param)
{
    .vcf_scan(file, ..., param=param) 
})

setMethod(scanVcf, c("TabixFile", "RangedData"),
    function(file, ..., param)
{
    .vcf_scan(file, ..., param=param)
})

setMethod(scanVcf, c("TabixFile", "GRanges"),
    function(file, ..., param)
{
    .vcf_scan(file, ..., param=param)
})

setMethod(scanVcf, c("TabixFile", "ScanVcfParam"),
    function(file, ..., param)
{
    result <- scanVcf(file, ..., fixed=vcfFixed(param),
                      info=vcfInfo(param), geno=vcfGeno(param),
                      samples=vcfSamples(param), param=vcfWhich(param))
    if (vcfTrimEmpty(param))
        lapply(result, function(rng) {
            rng[["GENO"]] <- Filter(Negate(is.null), rng[["GENO"]])
            rng
        })
    else 
        result
})

setMethod(scanVcf, c("TabixFile", "missing"),
    function(file, ..., param)
{
    callGeneric(file, ..., param=ScanVcfParam())
})

setMethod(scanVcf, c("character", "ScanVcfParam"),
    function(file, ..., param)
{
    ## no ranges
    if (0L == length(vcfWhich(param))) {
        .vcf_scan_character(file, ..., fixed=vcfFixed(param), 
            info=vcfInfo(param), geno=vcfGeno(param),
            samples=vcfSamples(param))
    } else {
    ## ranges
      callGeneric(TabixFile(file), ..., param=param)
    }
})

setMethod(scanVcf, c("character", "missing"),
    function(file, ..., param)
{
    callGeneric(file, ..., param=ScanVcfParam())
})

setMethod(scanVcf, c("connection", "missing"),
    function(file, ..., param)
{
    .vcf_scan_connection(file, ...)
})
