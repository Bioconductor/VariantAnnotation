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

.vcf_scan <- 
   function(file, ..., fixed=character(), info=character(),
            geno=character(), param)
{
    result <- tryCatch({
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }
        hdr <- scanVcfHeader(file)[[1]]
        samples <- hdr$Sample
        fmap <- .vcf_fixed(fixed)
        imap <- .vcf_map(hdr$Header[["INFO"]], info, nm="info")
        gmap <- .vcf_map(hdr$Header[["FORMAT"]], geno, nm="geno")
        tbxstate <- list(samples, fmap, imap, gmap)
        tbxsym <- getNativeSymbolInfo(".tabix_as_vcf",
                                      "VariantAnnotation")
        scanTabix(file, ..., param=param, tbxsym=tbxsym, tbxstate=tbxstate)
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             path(file), call. = FALSE)
    })

    .unpackVcf(result, hdr)
}

.vcf_scan_character <-
    function(file, ..., fixed=character(), info=character(),
             geno=character(), yieldSize=100000L)
{
    res <- tryCatch({
        file <- normalizePath(path.expand(file))
        if (!file.exists(file))
            stop("file does not exist")
        hdr <- scanVcfHeader(file)[[1]]
        samples <- hdr$Sample
        fmap <- .vcf_fixed(fixed)
        imap <- .vcf_map(hdr$Header[["INFO"]], info, "info")
        gmap <- .vcf_map(hdr$Header[["FORMAT"]], geno, "geno")
        result <- .Call(.scan_vcf_character, file,
                        as.integer(yieldSize), samples,
                        fmap, imap, gmap)
        names(result) <- "*:*-*"
        result
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             file, call. = FALSE)
    })

    .unpackVcf(res, hdr)
}

.vcf_scan_connection <-
    function(file, ..., fixed=character(), info=character(),
             geno=character())
{
    res <- 
        tryCatch({
            fl <- summary(file)$description
            hdr <- scanVcfHeader(fl)[[1]]
            samples <- hdr$Sample
            fmap <- .vcf_fixed(fixed)
            imap <- .vcf_map(hdr$Header[["INFO"]], info, "info")
            gmap <- .vcf_map(hdr$Header[["FORMAT"]], geno, "geno")

            txt <- readLines(file, ...)
            txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
            result <- .Call(.scan_vcf_connection, txt, samples,
                            fmap, imap, gmap)
            names(result) <- "*:*-*"
            result
        }, error = function(err) {
            stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
                 summary(file)$description, call. = FALSE)
        })

    .unpackVcf(res, hdr)
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
    result <- if (length(vcfWhich(param)) == 0) # no ranges
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
    if (length(vcfWhich(param)) == 0) {
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
                idx[idx == 0] <- d
                xrep <- rep(x, idx)
                x <- array(unlist(strsplit(xrep, ",", fixed=TRUE)),
                           dim=c(d, nrow(x), ncol(x)),
                           dimnames=list(NULL, NULL, colnames(x)))
                x <- aperm(x, c(2, 3, 1))
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
    if (length(x[[1]]$INFO) != 0) {
        info <- hdr[["Header"]][["INFO"]]
        if (is.null(info)) {
            warning("'INFO' vcf header info not found in file")
        } else {
            x <- lapply(x, function(elt, id, n, type) {
                elt[["INFO"]] <- .unpackVcfInfo(elt[["INFO"]], id, n, type)
                elt
            }, rownames(info), info$Number, info$Type)
        }
    }

    if (length(unlist(x[[1]]$GENO, use.names=FALSE)) != 0) {
        geno <- hdr[["Header"]][["FORMAT"]]
        if (is.null(geno)) {
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
