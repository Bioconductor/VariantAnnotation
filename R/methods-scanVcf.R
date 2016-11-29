### =========================================================================
### scanVcf methods 
### =========================================================================

selectSome <- S4Vectors:::selectSome
.vcf_usertag <-
    function(map, tag, nm, verbose)
{
    if (!identical(character(), tag)) {
        if (1L == length(tag) && is.na(tag)) {
            map[] <- NULL
        } else {
            ok <- tag %in% names(map)
            if (!all(ok)) {
                warning("ScanVcfParam ", sQuote(nm), " fields not found in ",
                        " header: ", paste(sQuote(tag[!ok]), collapse=" "))
                tag <- tag[ok]
            }
            map <- map[tag]  # user-requested order
        }
    }
    if (verbose) {
        msg <- paste0("found header lines for ", length(map), 
                      " ", sQuote(nm), " fields: ",
                      paste(selectSome(names(map)), collapse=", "))
        cat(msg, "\n")
    }
    if (!length(map) && nm == "info" && !length(tag))
        map <- list(list("1", character()))
    map
}

.vcf_map_fixed <-
    function(tag, verbose)
{
    map <- list(ALT=list("A", character()), QUAL=list("1", numeric()),
                FILTER=list("1", character()))
    c(list(rowRanges=NULL, REF=NULL), .vcf_usertag(map, tag, "fixed", verbose))
}

.vcf_map_samples <-
    function(fmt, tag, ...)
{
    map <- setNames(seq_along(fmt), fmt)
    if (identical(character(), tag)) {
        map                             # keep everyone
    } else if (isTRUE(is.na(tag))) {
        map[0]                          # drop everyone
    } else {
        if (any(nx <- !tag %in% fmt)) {
            warning("samples not in file: ",
                    paste(sQuote(tag[nx]), collapse=" "))
            tag <- tag[!nx]
        }
        map[] <- 0L                     # drop everyone, except...
        map[match(tag, fmt)] <- seq_along(tag) # these guys, and...
        map[rev(cumsum(rev(map)) > 0)]  # drop trailing unwanted samples
    }
}

## used for both INFO and GENO
.vcf_map <-
    function(fmt, tag, nm, verbose)
{
    numberOk <- grepl("^[AGR\\.[:digit:]+]$", fmt$Number)
    fmt$Number[!numberOk] <- "."
    mapType <- function(n, t) {
        if (t == "Flag") n <- "1"
        t <- switch(t, 
                    String=character(), 
                    Character=character(),
                    Integer=integer(), 
                    Float=numeric(), 
                    Flag=logical())
        list(n, t)
    }
    map <- Map(mapType, fmt$Number, fmt$Type)
    names(map) <- rownames(fmt)

    ## user selected
    .vcf_usertag(map, tag, nm, verbose)
}

.vcf_scan_header_maps <-
    function(file, fixed, info, geno, samples, verbose=FALSE)
{
    if (isTRUE(is.na(samples)))
       geno <- NA 
    hdr <- suppressWarnings(scanVcfHeader(file))
    fmap <- .vcf_map_fixed(fixed, verbose=verbose)
    imap <- .vcf_map(info(hdr), info, "info", verbose=verbose) 
    if (!is.null(names(imap))) {
        mapply(function(field, fname) {
                   if (field[[1]] == 0L && !is.logical(field[[2]]))
                       stop(paste0('only "flag" INFO fields should ',
                            'have Number = 0 in the header; ', fname, ' is "', 
                            class(field[[2]]), '"'))
                   else field 
               }, imap, names(imap)
        )
    }
    gmap <- .vcf_map(geno(hdr), geno, "geno", verbose)
    smap <- .vcf_map_samples(samples(hdr), samples)
    list(hdr=hdr, samples=smap, fmap=fmap, imap=imap, gmap=gmap)
}

.vcf_scan <- 
   function(file, ..., fixed=character(), info=character(),
            geno=character(), samples=character(), param,
            row.names=TRUE)
{
    tryCatch({
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }
        maps <- .vcf_scan_header_maps(file, fixed, info, geno, samples, ...)
        tbxstate <- maps[c("samples", "fmap", "imap", "gmap")]
        tbxsym <- getNativeSymbolInfo(".tabix_as_vcf", "VariantAnnotation")
        scanTabix(file, ..., param=param, tbxsym=tbxsym, 
                  tbxstate=tbxstate, row.names=row.names)
    }, scanTabix_io = function(err) {
            stop("scanVcf: ", conditionMessage(err), call. = FALSE)
    }, error = function(err) {
            stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
                 path(file), call. = FALSE)
    })
}

.vcf_scan_character <-
    function(file, ..., fixed=character(), info=character(),
             geno=character(), samples=character(), 
             yieldSize=100000L, row.names=TRUE)
{
    tryCatch({
        file <- normalizePath(path.expand(file))
        if (!file.exists(file))
            stop("file does not exist")
        maps <- .vcf_scan_header_maps(file, fixed, info, geno, samples, ...)
        result <- .Call(.scan_vcf_character, file, as.integer(yieldSize), 
                        maps$samples, maps$fmap, maps$imap, maps$gmap,
                        row.names)
        setNames(result, "*:*-*")
    }, scanTabix_io = function(err) {
            stop("scanVcf: ", conditionMessage(err), call. = FALSE)
    }, error = function(err) {
            stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
                 file, call. = FALSE)
    })
}

.vcf_scan_connection <-
    function(file, ..., fixed=character(), info=character(),
             geno=character(), samples=character(), row.names=TRUE)
{
    tryCatch({
        file <- summary(file)$description
        maps <- .vcf_scan_header_maps(file, fixed, info, geno, samples, ...)

        txt <- readLines(file, ...)
        txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
        result <- .Call(.scan_vcf_connection, txt, maps$samples,
                        maps$fmap, maps$imap, maps$gmap, row.names)
        setNames(result, "*:*-*")
    }, scanTabix_io = function(err) {
            stop("scanVcf: ", conditionMessage(err), call. = FALSE)
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
