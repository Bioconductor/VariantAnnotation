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
    rd <- rowData(obj)

    CHROM <- as.vector(seqnames(rd))
    POS <- start(rd)
    ID <- .makeVcfID(names(rd))
    REF <- as.character(ref(obj))
    ALT <- alt(obj)
    if (is(ALT, "DNAStringSetList")) {
        ALT <- as(ALT, "CharacterList")
    }
    ALT <- .pasteCollapse(ALT, ",")
    ALT[!nzchar(ALT)] <- "."
    QUAL <- qual(obj)
    QUAL[is.na(QUAL)] <- "."
    FILTER <- filt(obj)
    FILTER[is.na(FILTER)] <- "."
    INFO <- .makeVcfInfo(info(obj), length(rd))
    ans <- paste(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, sep = "\t")
    if (nrow(colData(obj)) > 0L) {
      GENO <- .makeVcfGeno(geno(obj))
      FORMAT <- GENO[,1]
      GENO <- GENO[,-1,drop=FALSE]
      genoPasted <- do.call(paste, c(split(GENO, col(GENO)), sep = "\t"))
      ans <- paste(ans, FORMAT, genoPasted, sep = "\t")
    }
    ans
}

.makeVcfID <- function(id, ...)
{
    if (is.null(id))
        "."
    else {
        idx <- grep(":", id, fixed=TRUE)
        id[idx] <- "."
        id
    }
}

.makeVcfFormatMatrix <- function(geno, cls, idx) {
  if (length(idx) > 0) {
    geno[idx] <- seqapply(geno[idx], function(x) {
      matrix(unlist(x, use.names=FALSE), nrow(x), prod(tail(dim(x), -1)))
    })
  }
 
  do.call(cbind, Map(function(elt, nms) {
### Should be discussed, but it seems like if we have a list matrix,
### we should look for elements that are empty, not a single NA.
### VO: If readVcf() was used, all empty fields were converted to NA.
    if (is.list(elt))
      haveData <- rowSums(matrix(elementLengths(elt), nrow(elt))) > 0
    else 
      haveData <- rowSums(!is.na(elt)) > 0
    as.character(ifelse(haveData, nms, NA_character_))
  }, as.list(geno), names(geno)))
}

.makeVcfFormat <- function(formatMat)
{
    if (ncol(formatMat) == 1L) {
        formatMat[is.na(formatMat)] <- ""
        .pasteCollapse(seqsplit(formatMat, row(formatMat)), ":")
    } else {
        keep <- !is.na(formatMat)
        .pasteCollapse(seqsplit(formatMat[keep], row(formatMat)[keep]), ":")
    }
}

.makeVcfGeno <- function(geno, ...)
{
    cls <- lapply(geno, class) 
    idx <- which(cls == "array")
    formatMat <- .makeVcfFormatMatrix(geno, cls, idx)
    FORMAT <- .makeVcfFormat(formatMat)

    nsub <- ncol(geno[[1]])
    nrec <- nrow(geno[[1]])
    arylst <- lapply(geno[idx], 
                  function(elt, nsub) 
                  {
                      m <- matrix(c(elt), ncol=nsub)
                      matrix(split(m, seq_len(nrec*nsub)), ncol=nsub)
                  }, nsub)
    geno[idx] <- SimpleList(arylst)

    genoMat <- matrix(unlist(as.list(geno), use.names = FALSE,
                             recursive = FALSE),
                      nsub * nrec, length(geno))

    ## convert NA values to '.' and get a simple character matrix
    genoMatFlat <- as.character(unlist(genoMat))
    genoMatFlat[is.na(genoMatFlat)] <- "."
    if (is.list(genoMat)) {
      genoMatList <- relist(genoMatFlat, PartitioningByEnd(genoMat))
      genoMatFlat <- .pasteCollapse(genoMatList, ",")
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
        keep <- !is.na(formatMatPerSub)
        genoListBySub <- seqsplit(genoMat[keep], row(genoMat)[keep])
        genoMatCollapsed <- matrix(.pasteCollapse(genoListBySub, ":"), nrec, nsub)
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
      collapsed <- .pasteCollapse(charList, ",")
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

    keep <- !is.na(infoMat)
    infoRows <- factor(row(infoMat), seq_len(nrow(infoMat)))
    infoList <- seqsplit(infoMat[keep], infoRows[keep])
    infoList[elementLengths(infoList) == 0L] <- "."
    .pasteCollapse(infoList, ";")
}

.pasteCollapse <- rtracklayer:::pasteCollapse
 
.makeVcfHeader <- function(obj, ...)
## FIXME : not writing out reference or sample
{
    hdr <- exptData(obj)[["header"]]
    header <- Map(.formatHeader, as.list(header(hdr)),
                  as.list(names(header(hdr))))
    samples <- samples(hdr)
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
            df$Description <- paste("\"", df$Description, "\"", sep="")
            prs <- paste(rep(colnames(df), each=nrow(df)), "=",
                         unlist(lapply(df, as.character), use.names=FALSE),
                         sep="")
            lst <- split(prs, seq_len(nrow(df)))
            lns <- .pasteCollapse(CharacterList(lst), collapse=",") 
            paste("##", nms, "=<ID=", rownames(df), ",", lns, ">", sep="")
        }
    }
}

