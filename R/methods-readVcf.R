### =========================================================================
### readVcf methods 
### =========================================================================

## TabixFile

setMethod(readVcf, c(file="TabixFile", param="ScanVcfParam"), 
    function(file, param, genome, ...)
{
    .readVcf(file, param, genome)
})

setMethod(readVcf, c(file="TabixFile", param="GRanges"),
    function(file, param, genome,  ...)
{
    .readVcf(file, param=ScanVcfParam(which=param), genome)
})

setMethod(readVcf, c(file="TabixFile", param="RangedData"),
    function(file, param, genome, ...)
{
    .readVcf(file, param=ScanVcfParam(which=param), genome)
})

setMethod(readVcf, c(file="TabixFile", param="RangesList"),
    function(file, param, genome, ...)
{
    .readVcf(file, param=ScanVcfParam(which=param), genome)
})

setMethod(readVcf, c(file="TabixFile", param="missing"), 
    function(file, param, genome, ...)
{
    .readVcf(file, genome=genome)
})

## character

setMethod(readVcf, c(file="character", param="ScanVcfParam"),
    function(file, param, genome, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, param, genome)
})

setMethod(readVcf, c(file="character", param="missing"),
    function(file, param, genome, ...)
{
    file <- .checkTabix(file)
    .readVcf(file, genome=genome)
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

.readVcf <- function(file, param, genome, ...)
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
                                  file, param=param, genome)
        } else {
            .scanVcfToVCF(scanVcf(file, param=param), file, genome)
        }
    }
}

## helpers

.scanVcfToLongGRanges <- function(vcf, file, param, genome, ...)
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

    ## fixed fields
    structural <- grep("<", vcf$ALT, fixed=TRUE)
    if (!is.null(vcf$ALT)) {
        if (!identical(integer(0), structural))
            ALT <- CharacterList(as.list(vcf$ALT))
        else
            ALT <- .toDNAStringSetList(vcf$ALT)
    } else {
        ALT <- vcf$ALT
    }
    if (!is.null(vcf$REF)) 
        REF <- .toDNAStringSet(vcf$REF)
    else
        REF <- vcf$REF

    if (!is.null(vcf$QUAL)) 
        QUAL <- DataFrame(vcf$QUAL)
    else
        QUAL <- vcf$QUAL
    if (!is.null(vcf$FILTER)) 
        FILTER <- DataFrame(vcf$FILTER)
    else
        FILTER <- vcf$FILTER
    fxfld <- list(REF=REF, ALT=ALT, QUAL=QUAL, FILTER=FILTER)
    fixed <- DataFrame(unlist(fxfld)) 
    dimnames(fixed) <- list(NULL, names(unlist(fxfld))) 

    ## rowData
    if (any(is.null(c(vcf$CHROM, vcf$POS, vcf$ID, vcf$REF)))) {
        rowData <- GRanges()
    } else if (!any(is.null(c(vcf$CHROM, vcf$POS, vcf$ID, vcf$REF)))) {
        rowData <- GRanges(seqnames=Rle(vcf$CHROM), ranges=IRanges(start=vcf$POS, 
            width=width(REF)))
        rowData <- .rowDataNames(vcf, rowData)
        genome(seqinfo(rowData)) <- genome
    }

    ## info 
    info <- .formatInfo(vcf$INFO,  hdr[["INFO"]])

    ## colData
    if (length(vcf$GENO) > 0) {
        samples <- colnames(vcf$GENO[[1]]) 
        colData <- DataFrame(Samples=seq_len(length(samples)), row.names=samples)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    VCF(rowData=rowData, colData=colData, exptData=SimpleList(HEADER=hdr), 
        fixed=fixed, info=info, geno=SimpleList(vcf$GENO))
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
