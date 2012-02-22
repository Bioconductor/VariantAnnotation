### =========================================================================
### readVcf methods 
### =========================================================================

## TabixFile

setMethod(readVcf, c(file="TabixFile", genome="character", 
          param="ScanVcfParam"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="GRanges"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="RangedData"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="RangesList"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character",
          param="missing"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome)
})

## character

setMethod(readVcf, c(file="character", genome="character",
          param="ScanVcfParam"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="character", genome="character",
          param="missing"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="character", genome="missing",
          param="missing"),
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing") 
})

.checkTabix <- function(x)
{
    if (1L != length(x)) 
        stop("'x' must be character(1)")
    if (grepl("\\.tbi$", x))
        TabixFile(sub("\\.tbi", "", x))
    else 
        x 
}

## .readVcf internal

.readVcf <- function(file, genome, param, ...)
{
    if (missing(param)) {
        .scanVcfToVCF(scanVcf(file), file, genome)
    } else {
        if (vcfAsGRanges(param)) {
            .scanVcfToLongGRanges(scanVcf(file, param=param),
                file, genome, param=param)
        } else {
            .scanVcfToVCF(scanVcf(file, param=param), file, genome, param)
        }
    }
}

## helpers

.scanVcfToLongGRanges <- function(vcf, file, genome, ...)
{
    if (1L < length(vcf)) {
        rnglen <- lapply(vcf, function(rng) length(rng[["rowData"]])) 
        rangeID <- as.factor(rep(names(vcf), rnglen))
        vcf <- lapply(names(vcf[[1]]), function(elt) {
            do.call(c, unname(lapply(vcf, "[[", elt)))
        })
        names(vcf) <- names(vcf[[1]]) 
    } else {
        rangeID <- as.factor(rep(names(vcf), length(vcf[[1]][["rowData"]])))
        vcf <- vcf[[1]]
    } 
    hdr <- scanVcfHeader(file)[[1]][["Header"]]

    ## single valued
    sdat <- list(REF=vcf[["REF"]], QUAL=vcf[["QUAL"]], FILTER=vcf[["FILTER"]])
    sdat <- DataFrame(rangeID=rangeID, ID=names(vcf[["rowData"]]), 
                      DataFrame(sdat[lapply(sdat, is.null) == FALSE])) 

    ## multi valued
    gdat <- .formatGeno(.combineLists(vcf[["GENO"]]))
    glen <- lapply(gdat, elementLengths) 
    nsmp <- ifelse(length(gdat) > 0, length(unique(gdat$SAMPLES)), 1)

    idat <- .formatInfo(.combineLists(vcf[["INFO"]]), hdr[["INFO"]])
    idat <- c(ALT=.formatALT(vcf[["ALT"]]), as.list(idat))
    ilen <- lapply(idat, function(x) rep(elementLengths(x), nsmp)) 
    idat <- lapply(idat, function(x) rep(x, nsmp)) 

    ## replicate field elements
    eltlen <- c(ilen, glen)
    if (identical(list(), eltlen)) {
        eltlen <- maxlen <- rep(1, nrow(sdat))
        multi <- list()
    } else {
        maxlen <- apply(do.call(cbind, eltlen), 1, prod)
        multi <- Map(function(elt, eltlen, maxlen) {
                     unlist(rep(elt, maxlen/eltlen), use.names=FALSE)
                 }, c(idat, gdat), eltlen, MoreArgs=list(maxlen))
    }
    single <- lapply(sdat, function(elt, maxlen, nsmp) {
                  rep(rep(elt, nsmp), maxlen)
              }, maxlen, nsmp)

    ## replicate rowData
    rowData <- vcf[["rowData"]]
    names(rowData) <- NULL
    genome(rowData) <- genome
    gr <- rep(rep(rowData, nsmp), maxlen)
    values(gr) <- DataFrame(c(single, multi)) 
    gr 
}

.scanVcfToVCF <- function(vcf, file, genome, ...)
{
    if (1L < length(vcf)) {
        rnglen <- lapply(vcf, function(rng) length(rng[["rowData"]]))
        rangeID <- as.factor(rep(names(vcf), rnglen))
        nms <- names(vcf[[1]])
        vcf <- lapply(nms, function(elt) {
            do.call(c, unname(lapply(vcf, "[[", elt)))
        })
        names(vcf) <- nms 
    } else {
        rangeID <- as.factor(rep(names(vcf), length(vcf[[1]][["rowData"]])))
        vcf <- vcf[[1]]
    }
    hdr <- scanVcfHeader(file)[[1]][["Header"]]

    ## rowData
    rowData <- vcf[["rowData"]]
    genome(seqinfo(rowData)) <- genome 
    values(rowData) <- DataFrame(rangeID=rangeID)


    ## fixed fields
    ALT <- .formatALT(vcf[["ALT"]])
    fx <- list(REF=vcf[["REF"]], ALT=ALT, QUAL=vcf[["QUAL"]],
               FILTER=vcf[["FILTER"]])
    fixed <- DataFrame(fx[lapply(fx, is.null) == FALSE]) 

    ## info 
    info <- .formatInfo(.combineLists(vcf[["INFO"]]), hdr[["INFO"]])

    ## colData
    if (length(vcf[["GENO"]]) > 0) {
        samples <- colnames(vcf[["GENO"]][[1]])
        colData <- DataFrame(Samples=seq_len(length(samples)),
                             row.names=samples)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    VCF(rowData=rowData, colData=colData, exptData=SimpleList(HEADER=hdr),
        fixed=fixed, info=info, geno=SimpleList(.combineLists(vcf[["GENO"]])))
}

.toDNAStringSet <- function(x)
{
    xx <- unlist(strsplit(x, ",", fixed = TRUE))
    ulist <- sub(".", "", xx, fixed = TRUE)
    ulist[is.na(ulist)] <- ""
    DNAStringSet(ulist)
}

.toDNAStringSetList <- function(x)
{
    dna <- .toDNAStringSet(x)
    idx <- elementLengths(strsplit(x, ",", fixed = TRUE))
    pbw <- PartitioningByWidth(idx)
    IRanges:::newCompressedList("DNAStringSetList",
        unlistData=dna, end=end(pbw))
}

.formatALT <- function(x)
{
    if (is.null(x))
        return(NULL)
    structural <- grep("<", x, fixed=TRUE)
    if (!identical(integer(0), structural))
        seqsplit(x, seq_len(length(x)))
    else
        .toDNAStringSetList(x)
}

.formatInfo <- function(x, hdr)
{
    if (length(x) == 0L)
        return(DataFrame())
    if (is.null(hdr)) {
        DF <- DataFrame(x)
        names(DF) <- names(x)
        return(DF)
    }
    ## matrices become compressed lists
    type <- hdr$Type[match(names(x), rownames(hdr))]
    idx <- which(lapply(x, function(elt) is.null(dim(elt))) == FALSE)
    if (0L != length(idx)) {
        for (i in idx) {
            x[[i]] <- .formatList(x[[i]], type[i])
        }
    }
    DF <- DataFrame(x)
    names(DF) <- names(x)
    DF
}

.formatList <- function(data, type)
{
    ncol <- ifelse(!is.na(dim(data)[3]), dim(data)[3], dim(data)[2])
    if (ncol > 1) 
       data <- split(unlist(data, use.names=FALSE), 
           rep(seq_len(dim(data)[1]), ncol))
 
    switch(type, 
        Integer = IntegerList(data),
        Float = NumericList(data),
        String = CharacterList(data),
        Logical = LogicalList(data))
}

.formatGeno <- function(x)
{
    if (length(x) == 0L)
        return(list())
    cls <- lapply(x, class)
    nvar <- dim(x[[1]])[1]
    nsmp <- dim(x[[1]])[2]
    nms <- colnames(x[[1]]) 

    ## collapse matrices and arrays
    for (i in which(cls == "array")) {
        dim(x[[i]]) <- c(nvar*nsmp, dim(x[[i]])[3])
        x[[i]] <- split(x[[i]], seq_len(nvar*nsmp))
    }
    for (i in which(cls == "matrix")) {
        dim(x[[i]]) <- c(nvar*nsmp, 1)
    }
    ## sample names become a data column
    c(list(SAMPLES=rep(nms, each=nvar)), x)
}

.combineLists <- function(x)
{
    ## combine elements across a list
    ## retain data structures
    sp <- split(x, unique(names(x)))
    lapply(sp, function(elt) {
        if (is(elt[[1]], "matrix")) {
            do.call(rbind, unname(elt))
        } else if (is(elt[[1]], "array")) { 
            d <- dim(elt[[1]])
            dat <- lapply(seq_len(d[3]), function(i)
                do.call(rbind, unname(lapply(elt, "[",,,i))))
            array(unlist(dat, use.names=FALSE),
                c(nrow(dat[[1]]), ncol(dat[[1]]), length(dat))) 
        } else {
            do.call(c, unname(elt))
        }
    })
}

