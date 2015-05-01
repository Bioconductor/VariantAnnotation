### =========================================================================
### Helper functions not exported
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readVcf()
###

.collapseLists <- function(vcf, param)
{
    idx <- sapply(vcf, function(elt) length(elt$rowRanges) > 0L)
    if (!sum(idx))
        return(vcf[[1]])
    if (length(vcf) > 1L)
        vcf <- vcf[idx]
    if (is(param, "ScanVcfParam"))
        paramRangeID <- names(unlist(vcfWhich(param), use.names=FALSE))[idx]
    else
        paramRangeID <- names(param)[idx]
    if (is.null(paramRangeID))
        paramRangeID <- rep(NA_character_, length(vcf))

    ## single range in 'which'
    if (1L == length(vcf)) {
        lst <- vcf[[1]]
        lst$paramRangeID <- as.factor(rep(paramRangeID, length(lst$rowRanges)))
    } else {
    ## multiple ranges in 'which'
        lst <- lapply(names(vcf[[1]]), function(elt) {
                   suppressWarnings(do.call(c, unname(lapply(vcf, "[[", elt))))
               })
        names(lst) <- names(vcf[[1]])
        len <- unlist(lapply(vcf, function(elt) length(elt$rowRanges)),
                      use.names=FALSE)
        paramRangeID <- as.factor(rep(paramRangeID, len))

        ## collapse info and geno
        info <- lst$INFO
        sp <- split(unname(info), unique(names(info)))
        sp <- sp[unique(names(info))]
        lst$INFO <-
            lapply(sp, function(elt) {
                d <- dim(elt[[1]])
                if (is(elt[[1]], "list")) {
                    as.matrix(elt)
                } else if (is(elt[[1]], "array") && !is.na(d[3])) {
                    pc <- lapply(seq_len(d[2]), function(i) {
                              do.call(rbind, lapply(elt, "[", ,i,))
                          })
                    array(do.call(c, pc),
                        c(length(lst$rowRanges), d[2], d[3]))
                } else {
                    do.call(c, elt)
                }
            })
        geno <- lst$GENO
        sp <- split(geno, unique(names(geno)))
        lst$GENO <-
            lapply(sp, function(elt) {
                d <- dim(elt[[1]])
                if (!is.na(d[3])) {
                    pc <- lapply(seq_len(d[3]), function(i) {
                              do.call(rbind, lapply(elt, "[", ,,i))
                          })
                    cmb <- array(do.call(c, pc),
                                 c(length(lst$rowRanges), d[2], d[3]))
                    cmb
                } else {
                    trans <- lapply(elt, t)
                    cmb <- matrix(do.call(c, trans), length(lst$rowRanges),
                                  d[2], byrow=TRUE)
                    cmb
                }
            })
        lst$paramRangeID <- paramRangeID
    }
    lst
}

.formatList <- function(data, type)
{
    switch(type,
        Integer = IntegerList(data),
        Float = NumericList(data),
        String = CharacterList(data),
        Character = CharacterList(data),
        Logical = LogicalList(data))
}

.toDNAStringSetList <- function(x)
{
    ### also used in predictCoding(), genotypeToSnpMatrix()
    pbw <- PartitioningByWidth(elementLengths(x))
    x <- unlist(x, use.names=FALSE)
    x[.isStructural(x)] <- ""
    xx <- sub(".", "", x, fixed=TRUE)
    relist(DNAStringSet(xx), pbw)
}

.formatALT <- function(x)
{
    if (is.null(x))
        return(NULL)
    if (any(.isStructural(unlist(x, use.names=FALSE))))
        CharacterList(x)
    else DNAStringSetList(x)
}

.isStructural <- function(x)
{
    grepl("<", x, fixed=TRUE) |
    grepl("[", x, fixed=TRUE) |
    grepl("]", x, fixed=TRUE)
}

.formatInfo <- function(x, hdr, nrecords)
{
    ## no data
    if (!length(x) || (1L == length(x) && all(x[[1]] == ".", na.rm=TRUE)))
        return(DataFrame(row.names=seq_len(nrecords)))
    ## data in file but not in header
    if (!length(hdr) && length(x[[1]])) {
        DF <- DataFrame(x)
        names(DF) <- names(x)
        return(DF)
    }
    ## matrices and arrays
    type <- hdr$Type[match(names(x), rownames(hdr))]
    idx <- which(lapply(x, is.array) == TRUE)
    if (0L != length(idx)) {
        for (i in idx) {
            dat <- x[[i]]
            d <- dim(dat)
            ncol <- ifelse(!is.na(d[3]), d[3], d[2])
            if (ncol > 1)
                 dat <- split(unlist(dat, use.names=FALSE),
                               rep(seq_len(d[1]), ncol))
            x[[i]] <- .formatList(dat, type[i])
        }
    }
    ## ragged lists
    lx <- which(lapply(x, is.list) == TRUE)
    if (0L != length(lx)) {
        for (i in lx) {
            x[[i]] <- .formatList(x[[i]], type[i])
        }
    }
    DF <- DataFrame(x)
    names(DF) <- names(x)
    DF
}

.rleRecycleVector <- function(x, len) {
  if (is(x, "Rle"))
    rep(x, length.out = len)
  else if (length(x) == 1L)
    Rle(x, len)
  else S4Vectors:::recycleVector(x, len)
}

.pasteCollapseRows <- function(x, sep = ",") {
  if (!is.matrix(x) || mode(x) != "character") {
    stop("'x' must be a matrix of mode character")
  }
  if (!isSingleString(sep) || nchar(sep) == 0L) {
    stop("'sep' must be a single, non-NA, non-empty string")
  }
  .Call(matrix_pasteCollapseRows, x, sep)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PolyPhen, SIFT and PROVEAN
###

.getKcol <- function(conn)
{
    sql <- "SELECT value FROM metadata WHERE name='Key column'"
    as.character(dbGetQuery(conn, sql))
}

.sqlIn <- function(vals)
{
    if (0L == length(vals))
        return("")
    sql <-
      lapply(seq_len(length(vals)),
        function(i) {
            v <- vals[[i]]
            if (!is.numeric(v))
            v <- paste("'", v, "'", sep="")
            v
        })
    paste(unlist(sql), collapse = ",")
}

.missingKeys <- function(x, keys, db)
{
    if (missing(keys))
        return(FALSE)
    if (any(mkeys <- !keys %in% keys(x))) {
        msg <- paste(BiocGenerics:::selectSome(keys[mkeys]), collapse=" ")
        warning(sum(mkeys), " keys not found in ", db, " database: ", msg,
                call.=FALSE)
    }
    all(mkeys)
}

.missingCols <- function(x, cols, db)
{
    if (missing(cols))
        return(FALSE)
    if (any(mcols <- !cols %in% columns(x))) {
        msg <- paste(BiocGenerics:::selectSome(cols[mcols], 5), collapse=" ")
        warning(sum(mcols), " columns not found in ", db, " database: ", msg,
                call.=FALSE)
        return(TRUE)
    } else {
        return(FALSE)
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isSNV helpers
###

.isSNV <- function(ref, alt) {
    nchar(ref) == 1L & nchar(alt) == 1L
}

.isDeletion <- function(ref, alt) {
    nchar(alt) == 1L & nchar(ref) > 1L & substring(ref, 1, 1) == alt
}

.isInsertion <- function(ref, alt) {
    nchar(ref) == 1L & 
    nchar(alt) > 1L & 
    substring(alt, 1, 1) == as.character(ref)
}

.isIndel <- function(ref, alt) {
    .isDeletion(ref, alt) | .isInsertion(ref, alt)
}

.isDelins <- function(ref, alt) {
    !.isIndel(ref, alt) & !.isSubstitution(ref, alt)
}

.isTransition <- function(ref, alt) {
    m1 <- match(as.character(ref), c("A", "G", "T", "C"))
    m2 <- match(as.character(alt), c("G", "A", "C", "T"))
    res <- (m1 - m2) == 0
    res[is.na(res)] <- FALSE
    res
}

.isSubstitution <- function(ref, alt) {
    nchar(ref) == nchar(alt)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### predictCoding()
###

.localCoordinates <- function(from, to, ignore.strand, ...)
{
    ## 'to' is a GRangesList of cds by transcript
    map <- mapToTranscripts(unname(from), to, ignore.strand=ignore.strand, ...)
    if (length(map) == 0) {
        res <- GRanges()
        mcols(res) <- DataFrame(REF=DNAStringSet(), ALT=DNAStringSetList(),
                                varAllele=DNAStringSet(), CDSLOC=IRanges(),
                                PROTEINLOC=IntegerList(), QUERYID=integer(),
                                TXID=character(), CDSID=IntegerList())
        return(res)
    }

    xHits <- map$xHits
    txHits <- map$transcriptsHits
    flat_to <- unlist(to) ## names needed for mapping

    ## FIXME: cdsid is expensive
    cdsid <- IntegerList(integer(0))
    map2 <- mapToTranscripts(unname(from)[xHits], flat_to,
                             ignore.strand=ignore.strand)
    cds <- mcols(flat_to)$cds_id[map2$transcriptsHits]
    if (length(cds)) {
        cdslst <- unique(splitAsList(cds, map2$xHits))
        cdsid <- cdslst
    }
    if (is.null(txid <- names(to)[txHits]))
        txid <- NA_integer_

    ## protein coordinates
    pends <- c(ceiling(start(map)/3), ceiling(end(map)/3))
    plocs <- unique(IntegerList(split(pends, rep(seq_len(length(pends)/2)), 2)))

    res <- from[xHits]
    strand(res) <- strand(map)
    mcols(res) <- append(values(res), 
        DataFrame(CDSLOC=ranges(map), 
                  PROTEINLOC=plocs, 
                  QUERYID=xHits, 
                  TXID=txid, CDSID=cdsid))
    res 
}

