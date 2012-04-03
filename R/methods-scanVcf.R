## scanVcf 

.vcf_usertag <-
    function(map, tag, nm, ...)
{
    if (!identical(character(), tag))
        if (1L == length(tag) && is.na(tag)) {
            map[] <- list(NULL)
        } else {
            map[!names(map) %in% tag] <- list(NULL)
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
    map <- list(ALT=character(), QUAL=numeric(), FILTER=character())
    c(list(rowData=NULL, REF=NULL), .vcf_usertag(map, tag, "fixed"))
}

.vcf_map <-
    function(fmt, tag, ...)
{
    ## map header FORMAT to R types
    map <- lapply(fmt$Type, switch,
                  String=character(), Integer=integer(),
                  Float=numeric(), Flag=logical())

    ## reassigns >1 and non-numeric
    d <- suppressWarnings(as.integer(fmt$Number))
    map[is.na(d) | d > 1] <- list(character())
    names(map) <- rownames(fmt)

    ## user selected
    .vcf_usertag(map, tag, ...)
}

.vcf_scan_header_maps <-
    function(file, fixed, info, geno)
{
    hdr <- scanVcfHeader(file)
    samples <- samples(hdr) 
    fmap <- .vcf_fixed(fixed)
    imap <- .vcf_map(info(hdr), info, nm="info")
    if (0L == length(imap))
        imap <- list(character())
    gmap <- .vcf_map(geno(hdr), geno, nm="geno")
    list(hdr=hdr, samples=samples, fmap=fmap, imap=imap, gmap=gmap)
}

.vcf_scan <- 
   function(file, ..., fixed=character(), info=character(),
            geno=character(), param)
{
    result <- tryCatch({
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }
        maps <- .vcf_scan_header_maps(file, fixed, info, geno)
        tbxstate <- maps[c("samples", "fmap", "imap", "gmap")]
        tbxsym <- getNativeSymbolInfo(".tabix_as_vcf",
                                      "VariantAnnotation")
        scanTabix(file, ..., param=param, tbxsym=tbxsym, tbxstate=tbxstate)
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             path(file), call. = FALSE)
    })

    .unpackVcf(result, maps$hdr)
}

.vcf_scan_character <-
    function(file, ..., fixed=character(), info=character(),
             geno=character(), yieldSize=100000L)
{
    res <- tryCatch({
        file <- normalizePath(path.expand(file))
        if (!file.exists(file))
            stop("file does not exist")
        maps <- .vcf_scan_header_maps(file, fixed, info, geno)
        result <- .Call(.scan_vcf_character, file,
                        as.integer(yieldSize), maps$samples,
                        maps$fmap, maps$imap, maps$gmap)
        names(result) <- "*:*-*"
        result
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             file, call. = FALSE)
    })

    .unpackVcf(res, maps$hdr)
}

.vcf_scan_connection <-
    function(file, ..., fixed=character(), info=character(),
             geno=character())
{
    res <- tryCatch({
        file <- summary(file)$description
        maps <- .vcf_scan_header_maps(file, fixed, info, geno)

        txt <- readLines(file, ...)
        txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
        result <- .Call(.scan_vcf_connection, txt, maps$samples,
                        maps$fmap, maps$imap, maps$gmap)
        names(result) <- "*:*-*"
        result
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             summary(file)$description, call. = FALSE)
    })

    .unpackVcf(res, maps$hdr)
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
    result <- if (0L == length(vcfWhich(param))) # no ranges
        scanVcf(path(file), ..., param=param)
    else                                # ranges
        scanVcf(file, ..., fixed=vcfFixed(param),
                info=vcfInfo(param), geno=vcfGeno(param),
                param=vcfWhich(param))
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
    callGeneric(path(file), ...)
})

setMethod(scanVcf, c("character", "ScanVcfParam"),
    function(file, ..., param)
{
    ## no ranges
    if (0L == length(vcfWhich(param))) {
        .vcf_scan_character(file, ..., fixed=vcfFixed(param), 
            info=vcfInfo(param), geno=vcfGeno(param))
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

## unpackVcf

.unpackVcfField <-
    function(x, id, n, type)
    ## convert sub-fields in 'x' to full R representation
{
    d <- suppressWarnings(as.integer(n))
    withCallingHandlers({
        x <- if (!is.na(d)) {
            if (1L < d) {
            ## >1 
                idx <- as.integer(!is.na(x))
                idx[0L == idx] <- d
                xrep <- rep(x, idx)
                x <- array(unlist(strsplit(xrep, ",", fixed=TRUE)),
                           dim=c(d, nrow(x), ncol(x)),
                           dimnames=list(NULL, NULL, colnames(x)))
                x <- aperm(sub("\\.", NA, x), c(2, 3, 1))
            }
            ## 1 and 0 
            switch(type, 
                   Flag=x,
                   Character=, String=x,
                   Integer={ mode(x) <- "integer"; x },
                   Float={ mode(x) <- "numeric"; x },
                   stop(sprintf("unhandled FORMAT type '%s'", type)))
        } else {
            ## non-numeric
            nrow <- nrow(x)
            dimnames <- dimnames(x)
            x <- strsplit(x, ",", fixed=TRUE)
            x <- switch(type, 
                        Character=, String=x,
                        Integer= lapply(x, as.integer),
                        Float=lapply(x, as.numeric),
                        stop(sprintf("unhandled FORMAT type '%s'", type)))
            matrix(x, nrow, dimnames=dimnames)
        }
    }, warning=function(w) {
        msg <- sprintf("unpackVcf field '%s': %s", id,
                       conditionMessage(w))
        warning(msg, call.=FALSE)
        invokeRestart("muffleWarning")
    })
}

.unpackVcfTag <- function(tag, id, n, type)
{
    if (is.null(names(tag)))
        stop("'tag' must be a named list")

    Map(function(elt, nm, id, n, type) {
        if (is.na(idx <- match(nm, id))) {
            msg <- sprintf("element '%s' not found in file header",
                           nm)
            stop(msg)
        }
        if (0L == length(elt) || is.null(elt))
            elt
        else
            .unpackVcfField(elt, id[idx], n[idx], type[idx])
    }, tag, names(tag), MoreArgs=list(id, n, type))
}

.unpackVcfInfo <-
    function(info, id, n, type)
{
    result <- .unpackVcfTag(info, id, n, type)
    lapply(result, function(elt) {
        if (is(elt, "list"))
            unlist(elt, recursive=FALSE, use.names=FALSE)
        else if (is(elt, "matrix") & ncol(elt) == 1)
            as.vector(elt)
        else 
            elt
    })
}

.unpackVcf <- function(x, hdr, ...)
{
    if (0L != length(x[[1]]$INFO)) {
        info <- info(hdr) 
        if (0L == length(info)) {
            warning("'INFO' vcf header info not found in file")
        } else {
            x <- lapply(x, function(elt, id, n, type) {
                elt[["INFO"]] <- .unpackVcfInfo(elt[["INFO"]], id, n, type)
                elt
            }, rownames(info), info$Number, info$Type)
        }
    }

    if (0L != length(unlist(x[[1]]$GENO, use.names=FALSE))) {
        geno <- geno(hdr) 
        if (0L == length(geno)) {
            warning("'FORMAT' vcf header info not found in file")
        } else {
            x <- lapply(x, function(elt, id, n, type) {
                elt[["GENO"]] <- .unpackVcfTag(elt[["GENO"]], id, n, type)
                elt
            }, rownames(geno), geno$Number, geno$Type)
        }
    }
    x
}
