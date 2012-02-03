### =========================================================================
### readVcf methods 
### =========================================================================

## TabixFile

setMethod(readVcf, c(file="TabixFile", genome="character", param="ScanVcfParam"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="GRanges"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="RangedData"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="RangesList"),
    function(file, genome, param, ...)
{
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="TabixFile", genome="character", param="missing"), 
    function(file, genome, param, ...)
{
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="TabixFile", genome="missing", param="ANY"), 
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing")
})


## character

setMethod(readVcf, c(file="character", genome="character", param="ScanVcfParam"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome, param)
})

setMethod(readVcf, c(file="character", genome="character", param="GRanges"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="character", genome="character", param="RangedData"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="character", genome="character", param="RangesList"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome, param=ScanVcfParam(which=param))
})

setMethod(readVcf, c(file="character", genome="character", param="missing"),
    function(file, genome, param, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome)
})

setMethod(readVcf, c(file="character", genome="missing", param="ANY"),
    function(file, genome, param, ...)
{
    stop("'genome' argument is missing")
})

.checkTabix <- function(x)
{
    if (grepl("[.]tbi$", x))
        TabixFile(gsub(".tbi", "", x, fixed=TRUE))
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
            p <- param
            if (identical(character(0), vcfInfo(param))) 
                slot(p, "info") <- NA_character_
            else if (identical(character(0), vcfGeno(param))) 
                slot(p, "geno") <- NA_character_
            .scanVcfToLongGRanges(scanVcf(file, param=p), file, genome, param=p)
        } else {
            .scanVcfToVCF(scanVcf(file, param=param), file, genome)
        }
    }
}

## helpers

.scanVcfToLongGRanges <- function(vcf, file, genome, param, ...)
{
    vcf <- vcf[[1]]
    ref <- .toDNAStringSet(vcf$REF)
    rowData <- GRanges(seqnames=Rle(vcf$CHROM), 
        ranges=IRanges(start=vcf$POS, width=width(ref)))
    rowData <- .rowDataNames(vcf, rowData)

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
                    elt[slen, ]
            } else {
                rep(elt, reps)
            }
        }, reps)

    ## FIXME: why can't DF go in list via lapply above
    idx <- which(lapply(lst, is.array) == TRUE) 
    if (length(idx) != 0)
        for (i in idx) 
            ll[[i]] <- DataFrame(I(matrix(ll[[i]], nrow=nrow(ll[[i]]))))
    DF <- DataFrame(ll)
    names(DF) <- names(lst)
    DF
}


.scanVcfToVCF <- function(vcf, file, genome, ...)
{
    vcf <- vcf[[1]]
    hdr <- scanVcfHeader(file)[[1]][["Header"]]

    ## fixed fields
    structural <- grep("<", vcf$ALT, fixed=TRUE)
    if (!identical(integer(0), structural))
        alt <- CharacterList(as.list(vcf$ALT))
    else
        alt <- .toDNAStringSetList(vcf$ALT)

    ref <- .toDNAStringSet(vcf$REF)
    fixedFields <- DataFrame(REF=ref, ALT=alt, QUAL=vcf$QUAL, FILTER=vcf$FILTER)

    ## rowData
    rowData <- GRanges(seqnames=Rle(vcf$CHROM), ranges=IRanges(start=vcf$POS, 
        width=width(ref)))
    rowData <- .rowDataNames(vcf, rowData)
    genome(seqinfo(rowData)) <- genome

    ## info 
    info <- .formatInfo(vcf$INFO,  hdr[["INFO"]])

    ## colData
    if (length(vcf$GENO) > 0) {
        samples <- colnames(vcf$GENO[[1]]) 
        colData <- DataFrame(Samples=seq_len(length(samples)), row.names=samples)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    ## FIXME : .unpack still returns lists for Number="."
    if (length(vcf$GENO) > 0)
        geno <- lapply(vcf$GENO, function(elt) {
                    if (is.list(elt))
                        do.call(rbind, elt)
                    else
                        elt
                    })
    else
        geno <- list()

    VCF(rowData=rowData, colData=colData, exptData=SimpleList(HEADER=hdr), 
        fixedFields=fixedFields, info=info, geno=SimpleList(geno))
}

.rowDataNames <- function(vcf, rowData, ...)
## create chrom:start rownames for records with no ID
{
    idx <- vcf$ID == "."
    vcf$ID[idx] <- paste("chr", seqnames(rowData[idx]), ":", 
        start(rowData[idx]), sep="") 
    names(rowData) <- vcf$ID
    rowData
}


.toDNAStringSet <- function(x)
{
    xx <- unlist(strsplit(x, ",", fixed = TRUE))
    ulist <- gsub(".", "", xx, fixed = TRUE)
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
    if (is.null(hdr))
        return(CharacterList(as.list(x[[1]])))

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
