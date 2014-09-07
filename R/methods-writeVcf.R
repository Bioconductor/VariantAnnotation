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
    FIXED <- paste(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, sep = "\t")
    if (!is.null(GENO <- .makeVcfGeno(geno(obj), dim(obj))))
        paste(FIXED, GENO, sep="\t")
    else
        FIXED 
}

.makeVcfGeno <- function(geno, dvcf, ...)
{
    if (dvcf[2] == 0L)
        return(NULL)
    if (any(is.null(names(geno)))) {
        warning("all 'geno(<VCF>)' fields must be named")
        return(NULL)
    }
    if ("GT" %in% names(geno)) {
        geno <- geno[c("GT", setdiff(names(geno), "GT"))]
    }
    .Call(.make_vcf_geno, names(geno), as.list(geno), c(":", ","), 
          dvcf, sapply(geno, function(x) dim(x)[3])) 
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
