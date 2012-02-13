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
            if (identical(character(0), vcfFixed(param))) 
                slot(param, "fixed") <- NA_character_
            else if (identical(character(0), vcfInfo(param))) 
                slot(param, "info") <- NA_character_
            else if (identical(character(0), vcfGeno(param))) 
                slot(param, "geno") <- NA_character_
            .scanVcfToLongGRanges(scanVcf(file, param=param),
                                  file, genome, param=param)
        } else {
            .scanVcfToVCF(scanVcf(file, param=param), file, genome)
        }
    }
}

## helpers

.scanVcfToLongGRanges <- function(vcf, file, genome, param, ...)
{
    vcf <- vcf[[1]]
    rowData <- vcf$rowData

    elts <- na.omit(c(vcfInfo(param), vcfGeno(param)))
    dat <- c(vcf$INFO[names(vcf$INFO) %in% elts],
         vcf$GENO[names(vcf$GENO) %in% elts])
    eltrep <- lapply(dat, function(elt, dim0) {
               if (is.list(elt))
                   elementLengths(elt)
               else
                   rep(1, dim0)
           }, dim0=length(rowData)) 
    maxrep <- apply(do.call(cbind, eltrep), 1, max)
    unwind <- .unwind(maxrep, dat)
    gr <- rowData[rep(seq_len(length(rowData)), maxrep)]
    values(gr) <- unwind 
    gr 
}

.unwind <- function(reps, lst, ...)
{
    ## replication of rows
    ll <- lapply(lst, 
        function(elt, reps) {
            if (is(elt, "list")) {
                mt <- elementLengths(elt) == reps 
                newrep <- rep(1, length(mt))
                newrep[mt == FALSE] <- reps[mt == FALSE]
                unlist(rep(elt, newrep), use.names=FALSE) 
            } else if (is(elt, "array")) {
                slen <- rep(seq_len(length(reps)), reps)
                if (length(dim(elt)) == 3)
                    matrix(elt, ncol=dim(elt)[3])[slen, ]
                else
                    matrix(elt[slen, ], ncol=ncol(elt), nrow=length(slen))
            } else {
                rep(elt, reps)
            }
        }, reps)

    ## FIXME: DF in list via lapply above
    idx <- which(lapply(lst, is.array) == TRUE) 
    if (length(idx) != 0)
        for (i in idx) 
            ll[[i]] <- DataFrame(I(matrix(ll[[i]], nrow=nrow(lst[[i]]))))
    DF <- DataFrame(ll)
    names(DF) <- names(lst)
    DF
}


.scanVcfToVCF <- function(vcf, file, genome, ...)
{
    vcf <- vcf[[1]]
    hdr <- scanVcfHeader(file)[[1]][["Header"]]

    ## rowData
    rowData <- vcf$rowData
    genome(seqinfo(rowData)) <- genome 

    ## fixed fields
    REF <- vcf$REF
    structural <- grep("<", vcf$ALT, fixed=TRUE)
    if (!identical(integer(0), structural))
        ALT <- seqsplit(vcf$ALT, seq_len(length(vcf$ALT)))
    else
        ALT <- .toDNAStringSetList(vcf$ALT)
    REF <- vcf$REF
    QUAL <- DataFrame(vcf$QUAL)
    FILTER <- DataFrame(vcf$FILTER)
    fixed <- DataFrame(REF, ALT, QUAL, FILTER)
    dimnames(fixed) <- list(NULL, c("REF", "ALT", "QUAL", "FILTER"))

    ## info 
    info <- .formatInfo(vcf$INFO, hdr[["INFO"]])

    ## colData
    if (length(vcf$GENO) > 0) {
        samples <- colnames(vcf$GENO[[1]]) 
        colData <- DataFrame(Samples=seq_len(length(samples)),
                             row.names=samples)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    VCF(rowData=rowData, colData=colData, exptData=SimpleList(HEADER=hdr), 
        fixed=fixed, info=info, geno=SimpleList(vcf$GENO))
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

.newCompressedList <- IRanges:::newCompressedList
.formatInfo <- function(x, hdr)
{
    if (length(x) == 0L)
        return(DataFrame())
    if (is.null(hdr)) {
        DF <- DataFrame(x)
        names(DF) <- names(x)
        return(DF)
    }
    type <- hdr$Type[match(names(x), rownames(hdr))]
    idx <- which(lapply(x, is.list) == TRUE)
    if (length(idx) != 0) {
        for (i in idx) {
            x[[i]] <- switch(type[[i]],
                Integer = .newCompressedList("CompressedIntegerList", x[[i]]),
                Float = .newCompressedList("CompressedNumericList", x[[i]]),
                String = .newCompressedList("CompressedCharacterList", x[[i]]),
                Logical = .newCompressedList("CompressedLogicalList", x[[i]]))
        }
    }
    idx <- which(lapply(x, is.array) == TRUE)
    if (length(idx) != 0)
        for (i in idx)
            x[[i]] <- DataFrame(I(matrix(x[[i]], ncol=dim(x[[i]])[3])))
    DF <- DataFrame(x)
    names(DF) <- names(x)
    DF
}

