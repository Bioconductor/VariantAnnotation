### =========================================================================
### vaFilter methods 
### =========================================================================

setMethod(.vaValidity, "VAFilter", function (object) {
    msg <- NULL
    fmls <- formals(object)
    if (length(fmls) == 0 || names(fmls)[[1]] != "x")
        msg <- c(msg, paste("'filter' must have at least one argument, 'x'"))
    if (is.null(msg)) TRUE else msg
})

setMethod(vaFilter, "missing", function(fun, name, ...) {
    vaFilter(function(x) !logical(length(x)), name=name, ...)
})

setMethod(vaFilter, "function",
          function(fun, name, ...) 
{
    name <- mkScalar(as.character(name))
    fmls <- formals(fun)
    if (length(fmls) == 0 || names(fmls)[[1]] != "x")
        stop("UserArgumentMismatch",
             "'filter' must have at least one argument, 'x'")

    env <- new.env(parent=environment(fun))
    env[[".stats"]] <- NULL
    fun <- eval(substitute(function(x, subset=TRUE) {
        res <- FUN(x)
        VAFilterResult(x, res, NAME, subset)
    }, list(FUN=fun, NAME=name)))
    environment(fun) <- env

    new("VAFilter", fun, name=name, ...)
})

setMethod(vaFilter, "VAFilter", function(fun, name, ...) {
    slot(fun, ".Data")
})

setMethod(name, "VAFilter", function(x, ...) slot(x, "name"))

setMethod(show, "VAFilter", function(object) {
    cat("class:", class(object), "\n")
    cat("name:", name(object), "\n")
    cat("use vaFilter(object) to see filter\n")
})

dbSNPFilter <-
    function(dbSNP=character(0), .name="dbSNPFilter")
{
    #.check_type_and_length(regex, "character", 0:1)
    require(dbSNP, quietly=TRUE, character.only=TRUE)
    vaFilter(function(x) {
        .idx <- logical(length(x))
        nms <- rep(runValue(seqnames(x)), runLength(seqnames(x)))
        df <- data.frame(chrom=nms, start=start(x), index=seq_len(length(x)))
        snps <- df[width(x) == 1,]
        chr <- names(getSNPcount())[names(getSNPcount()) %in% snps$chrom]
        res <- lapply(chr, function(i, snps) {
            s <- getSNPlocs(i)
            q <- snps[snps$chrom %in% i,]
            dbsnp <- q$start %in% s$loc
            q$index[dbsnp]
        }, snps)
        .idx[do.call(rbind, res)] <- TRUE
        .idx
    }, name = .name)
}

regionFilter <-
    function(txdb, region="coding", .name="regionFilter")
{
    #.check_type_and_length(min, "numeric", 1)
    vaFilter(function(x) {
        loc <- locateVariants(x, txdb)
        res <- logical(length(x))
        res[unique(loc$queryID[loc$location %in% region])] <- TRUE
        res
    }, name=.name)
}

compose <-
    function(filt, ..., .name)
{
    lst <- if (missing(filt)) list(...) else list(filt, ...)
    #for (`filt, ...` in lst)
    #    .check_type_and_length(`filt, ...`, "VAFilter", NA)
    if (missing(.name))
        .name <- paste(sapply(lst, name), collapse=" o ")
    vaFilter(function(x) {
        .idx <- VAFilterResult(x, !logical(length(x)), subset=FALSE)
        for (elt in rev(lst)){
            .idx <- .idx & elt(x, subset=FALSE)@.Data
        }
        .idx
    }, name = .name)
}

