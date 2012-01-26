## scanVcf 

.vcf_map <-
    function(fmt, tag, ...)
{
    ## map header FORMAT to R types
    map <- lapply(fmt$Type, switch,
                  String=character(), Integer=integer(),
                  Float=numeric(), Flag=logical())
    ## FIXME: should these be parsed to true data type in C 
    ##        (reassigns >1 and non-numeric))
    d <- suppressWarnings(as.integer(fmt$Number))
    map[is.na(d)] <- list(character())
    map[d > 1] <- list(character())
    names(map) <- rownames(fmt)

    ## user selected 
    if (!identical(character(), tag))
        if (1L == length(tag) && is.na(tag)) {
            map[] <- list(NULL)
        } else {
            map[!names(map) %in% tag] <- list(NULL)
            if (!all(tag %in% names(map)))
                warning(paste("values in ScanVcfParam(", tag, "=c(...)) ",
                        "not present in header :", 
                        tag[!tag %in% names(map)], sep=""))
        }
    map
}

.vcf_scan <-
    function(file, ..., info=character(), geno=character(), param)
{
    res <- 
        tryCatch({
            if (!isOpen(file)) {
                open(file)
                on.exit(close(file))
            }
            hdr <- scanVcfHeader(file)[[1]]
            samples <- hdr$Sample
            imap <- .vcf_map(hdr$Header[["INFO"]], info)
            gmap <- .vcf_map(hdr$Header[["FORMAT"]], geno)
            tbx <- scanTabix(file, param=param) 
            result <- .Call(.scan_vcf, tbx, samples, imap, gmap)
            #names(result) <- sprintf("%s:%d-%d", space, start, end)
            result
        }, error = function(err) {
            stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
                 path(file), call. = FALSE)
        })

    .unpackVcf(res, hdr)
}

.vcf_scan_connection <-
    function(file, ..., info=character(), geno=character())
{
    res <- 
        tryCatch({
            fl <- summary(file)$description
            hdr <- scanVcfHeader(fl)[[1]]
            samples <- hdr$Sample
            imap <- .vcf_map(hdr$Header[["INFO"]], info)
            gmap <- .vcf_map(hdr$Header[["FORMAT"]], geno)

            txt <- readLines(file, ...)
            txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
            result <- .Call(.scan_vcf_connection, txt, samples, imap, gmap)
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
    ## no ranges
    if (length(vcfWhich(param)) == 0) 
        result <- callGeneric(path(file), param=param)
    else 
    ## ranges
        result <- scanVcf(file, ..., info=vcfInfo(param), geno=vcfGeno(param), 
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
    callGeneric(path(file))
})

setMethod(scanVcf, c("character", "ScanVcfParam"),
    function(file, ..., param)
{
    ## no ranges
    if (length(vcfWhich(param)) == 0) {
        con <- file(file)
        on.exit(close(con))
        .vcf_scan_connection(con, info=vcfInfo(param), geno=vcfGeno(param))
    } else {
    ## ranges
      callGeneric(TabixFile(file), ..., param=param)
    }
})

setMethod(scanVcf, c("character", "missing"),
    function(file, ..., param)
{
    con <- file(file)
    on.exit(close(con))
    callGeneric(con, ...)
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
            x <- apply(x, 1, function(i) strsplit(i, ",", fixed=TRUE)) 
            x <- switch(type, 
                        Character=, String=x,
                        Integer=lapply(x, lapply, as.integer),
                        Float=lapply(x, lapply, as.numeric),
                        stop(sprintf("unhandled FORMAT type '%s'", type)))
            names(x) <- NULL
            x
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
    ## handle NULL elements
        if (is.null(elt))
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

.unpackVcfGeno <-
    function(geno, id, n, type)
{
    .unpackVcfTag(geno, id, n, type)
    #lapply(result, function(elt) {
    #    if (is(elt, "list"))
    #        unlist(elt, recursive=FALSE, use.names=FALSE)
    #    else 
    #        elt
    #})
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
                elt[["GENO"]] <- .unpackVcfGeno(elt[["GENO"]], id, n, type)
                elt
            }, rownames(geno), geno$Number, geno$Type)
        }
    }
    x
}
