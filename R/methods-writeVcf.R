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
    GENO <- .makeVcfGeno(geno(obj))
    mat <- paste(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, GENO, sep = "\t")
    stopifnot(length(mat) == length(rd))
    mat
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

## FORMAT names listed once per row, specific to variant not sample
.makeVcfFormat <- function(geno)
{
    flst <- Map(function(g, nms) {
                    d <- dim(g)
                    sdim <- ifelse (is.na(d[3]), d[2], d[2]*d[3])
                    tapply(as.vector(g), rep(seq_len(d[1]), sdim), 
                        function(x)
                            ifelse (all(is.na(x)), "", nms))
                }, g=geno, nms=as.list(names(geno)))
    fvec <- do.call(paste, c(flst, sep=":"))
    s1 <- gsub("(::)", ":", fvec)
    s2 <- gsub("(^:|:$)", "", s1)
    s2
}

.makeVcfGeno <- function(geno, ...)
{
    if (length(geno) == 0L)
        return(NULL)

    fmt <- .makeVcfFormat(geno)
    ## collapse variable dimensions
    dlst <- lapply(geno, 
        function(g) {
            if (!is.na(dim(g)[3])) {
                tmp <- as.character(g[,,1]) 
                for (i in seq_len(dim(g)[3] - 1) + 1)
                    tmp <- paste(tmp, g[,,i], sep=",")
            } else {
                tmp <- as.character(g)
            }
        tmp 
        }
    )
 
    ## collapse across variables
    dvec <- do.call(paste, c(dlst, sep=":"))
    dcmb <- c(fmt, gsub("NA", "", dvec, fixed=TRUE))
    s1 <- gsub("(::)", ":", dcmb)
    s2 <- gsub("(^:|:$)", "", s1)
    dsep <- split(s2, rep(seq_along(fmt), ncol(geno[[1]]) + 1))
    do.call(rbind, lapply(dsep, function(x) paste(x, collapse="\t")))
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

