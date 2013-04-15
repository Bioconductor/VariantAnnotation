### =========================================================================
### Helper functions not exported 
### =========================================================================

.collapseLists <- function(vcf, param)
{
    idx <- sapply(vcf, function(elt) length(elt$rowData) > 0L)
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
        vcf[[1]][["paramRangeID"]] <- 
            as.factor(rep(paramRangeID, length(vcf[[1]][["rowData"]]))) 
        vcf[[1]]
    } else {
    ## multiple ranges in 'which'
        lst <- lapply(names(vcf[[1]]), function(elt) {
                   suppressWarnings(do.call(c, unname(lapply(vcf, "[[", elt))))
               })
        names(lst) <- names(vcf[[1]])
        len <- unlist(lapply(vcf, function(elt) length(elt$rowData)),
                      use.names=FALSE)
        paramRangeID <- as.factor(rep(paramRangeID, len)) 

        ## collapse info and geno
        info <- lst$INFO
        sp <- split(unname(info), unique(names(info)))
        sp <- sp[unique(names(info))] 
        lst[["INFO"]] <- 
            lapply(sp, function(elt) {
                d <- dim(elt[[1]])
                if (is(elt[[1]], "list")) {
                    as.matrix(elt)
                } else if (is(elt[[1]], "array") && !is.na(d[3])) {
                    pc <- lapply(seq_len(d[2]), function(i) {
                              do.call(rbind, lapply(elt, "[", ,i,)) 
                          })
                    array(do.call(c, pc), 
                        c(length(lst[["rowData"]]), d[2], d[3]))
                } else {
                    do.call(c, elt)
                }
            }) 
        geno <- lst[["GENO"]]
        sp <- split(geno, unique(names(geno)))
        lst[["GENO"]] <- 
            lapply(sp, function(elt) {
                d <- dim(elt[[1]])
                if (!is.na(d[3])) {
                    pc <- lapply(seq_len(d[3]), function(i) {
                              do.call(rbind, lapply(elt, "[", ,,i)) 
                          })
                    cmb <- array(do.call(c, pc), 
                                 c(length(lst[["rowData"]]), d[2], d[3]))
                    cmb
                } else {
                    trans <- lapply(elt, t)
                    cmb <- matrix(do.call(c, trans), length(lst[["rowData"]]), 
                                  d[2], byrow=TRUE)
                    cmb
                }
            }) 
        lst$paramRangeID <- paramRangeID
        lst
    }
}

.formatList <- function(data, type)
{
    switch(type, 
        Integer = IntegerList(data),
        Float = NumericList(data),
        String = CharacterList(data),
        Logical = LogicalList(data))
}

.toDNAStringSetList <- function(x)
{
    if (is(x, "CharacterList")) {
        ## Called from predictCoding, genotypeToSnpMatrix, etc.
        pbw <- PartitioningByWidth(elementLengths(x))
        x <- unlist(x, use.names=FALSE)
        str <- .isStructural(x)
        x[str] <- "" 
    } else {
        ## Called from .formatALT
        x <- strsplit(x, ",", fixed=TRUE)
        pbw <- PartitioningByWidth(elementLengths(x))
        x <- unlist(x, use.names=FALSE)
    }
    xx <- sub(".", "", x, fixed=TRUE)
    relist(DNAStringSet(xx), pbw)
}

.formatALT <- function(x)
{
    if (is.null(x))
        return(NULL)
    str <- .isStructural(x) 
    if (sum(str) != 0L) {
        lst <- strsplit(unlist(x, use.names=FALSE), ",")
        lst[str] <- x[str] 
        new("CompressedCharacterList",
            unlistData=unlist(lst, use.names=FALSE),
            partitioning=PartitioningByEnd(lst))
    } else {
        .toDNAStringSetList(x)
    }
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
    if (0L == length(x) ||
       (1L == length(x) && all(names(x) == "INFO")))
        return(DataFrame(row.names=seq_len(nrecords)))
    ## data in file but not in header
    if (0L == length(hdr)) {
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


.listCumsum <- function(x) {
  x_unlisted <- unlist(x, use.names=FALSE)
  x_cumsum <- cumsum(x_unlisted)
  x_part <- PartitioningByWidth(elementLengths(x))
  x_cumsum - rep(x_cumsum[start(x_part)] - x_unlisted[start(x_part)],
                 width(x_part))
}

.listCumsumShifted <- function(x) {
  cs <- .listCumsum(x)
  shifted <- c(0L, head(cs, -1))
  shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
  shifted
}

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

.setActiveSubjectSeq <-
    function(query, subject)
    ## set active state of circular sequences of subject to FALSE;
    ## warn if query contains circular sequences
{
    queryseq <- seqlevels(query)
    circular <- isCircular(subject)
    circNames <- intersect(queryseq, names(circular)[circular])
    if (0L != length(circNames))
        warning("circular sequence(s) in query '",
                paste(circNames, sep="' '"), "' ignored")

    isActiveSeq(subject)[] <- FALSE
    isActiveSeq(subject)[setdiff(queryseq, circNames)] <- TRUE
}
