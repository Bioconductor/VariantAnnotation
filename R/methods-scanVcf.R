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
            if (!all(tag %in% names(map))) {
                warning("ScanVcfParam '", nm, "' fields not present: '",
                        paste(tag[!tag %in% names(map)], collapse="' ' "),
                        "'")
            } else {
                fac <- match(c(tag, names(map)[!names(map) %in% tag]), names(map))
                map <- map[fac]
            }
        }
    map
}

.vcf_map_fixed <-
    function(tag, ...)
{
    map <- list(ALT=list("A", character()), QUAL=list("1", numeric()),
                FILTER=list("1", character()))
    c(list(rowData=NULL, REF=NULL), .vcf_usertag(map, tag, "fixed"))
}

.vcf_map_samples <-
    function(fmt, tag, ...)
{
    if (identical(character(), tag)) { 
        smap <- seq_along(fmt)
    } else if (isTRUE(is.na(tag))) {
        smap <- rep(0L, length(fmt))
    } else {
        if (any(nx <- !tag %in% fmt)) {
            warning(paste("samples", sQuote(tag[nx]), "not found in file"))
            tag <- tag[!nx]
        }
        smap <- rep(0L, length(fmt))
        smap[match(tag, fmt)] <- seq_along(tag)
    } 
    names(smap) <- fmt
    smap
}

.vcf_map <-
    function(fmt, tag, ...)
{
    numberOk <- grepl("^[AG\\.[:digit:]+]$", fmt$Number)
    fmt$Number[!numberOk] <- "."
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
    if (isTRUE(is.na(samples)))
       geno <- NA 
    hdr <- suppressWarnings(scanVcfHeader(file))
    smap <- .vcf_map_samples(samples(hdr), samples)
    gmap <- .vcf_map(geno(hdr), geno, nm="geno")
    fmap <- .vcf_map_fixed(fixed)
    imap <- .vcf_map(info(hdr), info, nm="info")
    if (0L == length(imap))
        imap <- list(list("1", character()))
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
