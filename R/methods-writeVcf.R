### =========================================================================
### writeVcf methods
### =========================================================================

setMethod(writeVcf, c("VCF", "character"),
    function(obj, filename, index = FALSE, ...)
{
    con <- file(filename, open="w")
    on.exit(close(con))

    writeVcf(obj, con, index=index, ...)
})

setMethod(writeVcf, c("VCF", "connection"),
    function(obj, filename, index = FALSE, ...)
{
    if (!isTRUEorFALSE(index))
        stop("'index' must be TRUE or FALSE")

    if (!isOpen(filename)) {
        open(filename)
        on.exit(close(filename))
    }

    scon <- summary(filename)
    headerNeeded <- !((scon$mode == "a") &&
                     file.exists(scon$description) &&
                     (file.info(scon$description)$size !=0))
    if (headerNeeded) {
        hdr <- .makeVcfHeader(obj)
        writeLines(hdr, filename)
    }

    if (index)
        obj <- sort(obj)
    mat <- .makeVcfMatrix(obj)
    writeLines(mat, filename)
    flush(filename)

    if (index) {
        filenameGZ <- bgzip(scon$description, overwrite = TRUE)
        indexTabix(filenameGZ, format = "vcf")
        unlink(scon$description)
        invisible(filenameGZ)
    } else {
        invisible(scon$description)
    }
})

.makeVcfMatrix <- function(obj)
{
    ## empty
    if (length(rd <- rowData(obj)) == 0)
        return(character())

    CHROM <- as.vector(seqnames(rd))
    POS <- start(rd)
    if (is.null(ID <- names(rd)))
        ID <- "."
    REF <- as.character(ref(obj))
    if (is.null(ALT <- alt(obj)))
        ALT <- rep(".", length(REF))
    if (is(ALT, "XStringSetList")) {
        ALT <- as(ALT, "CharacterList")
    }
    ALT <- as.character(unstrsplit(ALT, ","))
    ALT[nchar(ALT) == 0L | is.na(ALT)] <- "."
    if (is.null(QUAL <- qual(obj)))
        QUAL <- "."
    else
        QUAL[is.na(QUAL)] <- "."
    if (is.null(FILTER <- filt(obj)))
        FILTER <- "."
    else
        FILTER[is.na(FILTER)] <- "."
    INFO <- .makeVcfInfo(info(obj), length(rd))
    ans <- paste(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, sep = "\t")
    if (nrow(colData(obj)) > 0L) {
      GENO <- .makeVcfGeno(geno(obj))
      if (!is.null(GENO)) {
          FORMAT <- GENO[,1]
          GENO <- GENO[,-1,drop=FALSE]
          genoPasted <- do.call(paste, c(split(GENO, col(GENO)), sep = "\t"))
          ans <- paste(ans, FORMAT, genoPasted, sep = "\t")
      }
    }
    ans
}

.makeVcfFormatMatrix <- function(geno, cls, idx)
{
    if (length(geno) == 1L)
        return(matrix(rep(names(geno), nrow(geno[[1]]))))

    if (length(idx) > 0) {
        geno[idx] <- as(lapply(geno[idx],
            function(x) {
                matrix(unlist(x, use.names=FALSE), nrow(x), prod(tail(dim(x), -1)))
            }), "List")
    }

    do.call(cbind, Map(function(elt, nms) {
### Should be discussed, but it seems like if we have a list matrix,
### we should look for elements that are empty, not a single NA.
### VO: If readVcf() was used, all empty fields were converted to NA.
        if (is.list(elt))
          haveData <- rowSums(matrix(elementLengths(elt), nrow(elt))) > 0
        else
          haveData <- rowSums(!is.na(elt)) > 0
        ifelse(haveData, as.character(nms), NA_character_)
        }, as.list(geno), names(geno))
    )
}

.makeVcfFormat <- function(formatMat)
{
    .pasteCollapseRows(formatMat, ":")
}

.makeVcfGeno <- function(geno, ...)
{
    cls <- lapply(geno, class)
    idx <- which(cls == "array")
    formatMat <- .makeVcfFormatMatrix(geno, cls, idx)
    FORMAT <- .makeVcfFormat(formatMat)
    if (sum(nchar(FORMAT)) == 0L) {
        warning("all geno(<VCF>) fields are NA")
        return(NULL)
    }

    nsub <- ncol(geno[[1]])
    nrec <- nrow(geno[[1]])
    arylst <- lapply(geno[idx],
                  function(elt, nsub)
                  {
### FIXME: a bit inefficient, because we will generate elements like ".,."
###        for rows that are all NA, only to exclude them later.
                      elt[is.na(elt)] <- "."
                      apply(elt, 2, .pasteCollapseRows, ",")
                  }, nsub)
    geno[idx] <- SimpleList(arylst)

    genoMat <- matrix(unlist(as.list(geno), use.names = FALSE,
                             recursive = FALSE),
                      nsub * nrec, length(geno))

    ## handle Rle
    rle <- sapply(geno, function(x) is(x[[1]], "Rle"))
    if (any(rle))
        genoMat[,rle] <- lapply(genoMat[,rle], as.vector)

    ## convert NA values to '.'
    genoMatFlat <- as.character(unlist(genoMat, use.names=FALSE))
    genoMatFlat[is.na(genoMatFlat)] <- "."
    if (is.list(genoMat)) {
      genoMatList <- relist(genoMatFlat, PartitioningByEnd(genoMat))
      genoMatFlat <- unstrsplit(genoMatList, ",")
    }
    genoMat <- matrix(genoMatFlat, nrow(genoMat), ncol(genoMat))

    if (ncol(genoMat) == 1L) {
        if (ncol(genoMat) != length(FORMAT))
            cbind(FORMAT, matrix(genoMat, nrec, nsub))
        else
            cbind(FORMAT, genoMat)
    } else {
        formatMatPerSub <- matrix(rep(t(formatMat), nsub), nsub*nrec,
                                  length(geno), byrow=TRUE)
        genoMat[is.na(formatMatPerSub)] <- NA_character_
        genoMatCollapsed <- matrix(.pasteCollapseRows(genoMat, ":"), nrec, nsub)
        cbind(FORMAT, genoMatCollapsed)
    }

}

.makeVcfInfo <- function(info, nrecords, ...)
{
    if (ncol(info) == 0) {
      return(rep.int(".", nrecords))
    }

    ## Replace NA with '.' in columns with data.
    ## Columns with no data are set to NA.
    lists <- sapply(info, function(elt)
        is.list(elt) || is(elt, "List"))
    info[lists] <- lapply(info[lists], function(l) {
      charList <- as(l, "CharacterList")
      charList@unlistData[is.na(charList@unlistData)] <- "."
      collapsed <- unstrsplit(charList, ",")
      ifelse(sum(!is.na(l)) > 0L, collapsed, NA_character_)
    })

    ## Add names to non-NA data.
    infoMat <- matrix(".", nrow(info), ncol(info))
    logicals <- sapply(info, is.logical)
    infoMat[,logicals] <- unlist(Map(function(l, nm) {
      ifelse(l, nm, NA_character_)
    }, info[logicals], as(names(info)[logicals], "List")))

    infoMat[,!logicals] <- unlist(Map(function(i, nm) {
      ifelse(!is.na(i), paste0(nm, "=", i), NA_character_)
    }, info[!logicals], as(names(info)[!logicals], "List")))

    infoVector <- .pasteCollapseRows(infoMat, ";")
    infoVector[!nzchar(infoVector)] <- "."
    infoVector
}

.makeVcfHeader <- function(obj, ...)
{
    hdr <- header(obj)
    header <- Map(.formatHeader, as.list(header(hdr)),
                  as.list(names(header(hdr))))
    header <- c(header, .contigsFromSeqinfo(seqinfo(obj)))
    samples <- colnames(obj)
    colnms <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    if (length(geno(obj)) > 0L) {
      colnms <- c(colnms, "FORMAT", samples[!is.null(samples)])
    }
    colnms <- paste(colnms, collapse="\t")
    unlist(c(header, colnms), use.names=FALSE)
}

.formatHeader <- function(df, nms)
{
    if (nms == "META") {
        fd <- format(Sys.time(), "%Y%m%d")
        if ("fileDate" %in% rownames(df))
            df[rownames(df) == "fileDate", ] <- fd
        else
            df <- rbind(df, DataFrame(Value=fd, row.names="fileDate"))
        paste("##", rownames(df), "=", df[,1], sep="")
    } else {
        if ("Description" %in% colnames(df)) {
            if (nrow(df) == 0L)
                return(character())
            df$Description <-
              ifelse(is.na(df$Description), "\".\"",
                           paste("\"", df$Description, "\"", sep=""))
            df <- DataFrame(ID = rownames(df), df)
            prs <- paste(rep(colnames(df), each=nrow(df)), "=",
                         unlist(lapply(df, as.character), use.names=FALSE),
                         sep="")
            lst <- split(prs, row(df))
            lns <- unstrsplit(lst, ",")
            paste("##", nms, "=<", lns, ">", sep="")
        }
    }
}

.contigsFromSeqinfo <- function(si) {
  contig <- paste0("##contig=<ID=", seqnames(si))
  contig[!is.na(seqlengths(si))] <-
    paste0(contig, ",length=", seqlengths(si))[!is.na(seqlengths(si))]
  contig[!is.na(genome(si))] <-
    paste0(contig, ",assembly=\"", genome(si), "\"")[!is.na(genome(si))]
  paste0(contig, ">")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### VRanges methods
###

setMethod(writeVcf, "VRanges",
          function(obj, filename, ...)
          {
            writeVcf(as(obj, "VCF"), filename, ...)
          })
