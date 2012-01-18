### =========================================================================
### Helper functions not exported 
### =========================================================================

## for readVcf :

.scanVcfToLongGRanges <- function(vcf, file, genome, param, ...)
{
    vcf <- vcf[[1]]
    if (vcfAsGRanges(param) == "info")
        dat <- vcf$INFO
    else if (vcfAsGRanges(param) == "geno")
        dat <- vcf$GENO
    else
        stop("asGRanges must be specified as 'info' or 'geno'")
 
    ref <- .toDNAStringSet(vcf$REF)
    rowData <- GRanges(seqnames=Rle(vcf$CHROM), 
        ranges=IRanges(start=vcf$POS, width=width(ref)))
    rowData <- .rowDataNames(vcf, rowData)

    eltrep <- lapply(dat, function(elt, dim0) {
               if (is.list(elt))
                   elementLengths(elt)
               else
                   rep(1, dim0)
           }, dim0=length(rowData)) 
    maxrep <- apply(do.call(cbind, eltrep), 1, max)
    unwind <- .unwind(length(rowData), maxrep, dat)
    gr <- rowData[rep(seq_len(length(rowData)), maxrep)]
    values(gr) <- unwind 
    gr 
}

.unwind <- function(rep0, rep1, lst, ...)
{
   ll <- lapply(lst, function(elt, rep1) {
              if (is(elt, "list")) {
                  mt <- elementLengths(elt) == rep1 
                  newrep <- rep(1, length(mt))
                  newrep[mt == FALSE] <- rep1[mt == FALSE]
                  unlist(rep(elt, newrep), use.names=FALSE) 
              } else if (is(elt, "array")) {
                  slen <- rep(seq_len(length(rep1)), rep1)
                  matrix(elt, ncol=dim(elt)[3])[slen, ]
              } else {
                  rep(elt, rep1)
              }
          }, rep1)

    ## FIXME: why can't DF go in list via lapply above
    idx <- which(lapply(lst, class) == "array") 
    if (length(idx) != 0)
        for (i in idx) 
            ll[[i]] <- DataFrame(I(ll[[i]]))
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

    ## geno
    if (length(vcf$GENO) > 0) 
        geno <- lapply(vcf$GENO, function(elt) {
                    if (is.list(elt))
                        do.call(rbind, elt) 
                    else
                        elt
                    })
    else 
        geno <- list() 

    ## info 
    info <- .formatInfo(vcf$INFO,  hdr[["INFO"]])

    ## colData
    if (length(vcf$GENO) > 0) {
        samples <- colnames(geno[[1]]) 
        colData <- DataFrame(Samples=seq_len(length(samples)), row.names=samples)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    VCF(rowData=rowData, colData=colData, exptData=SimpleList(HEADER=hdr), 
        fixedFields=fixedFields, info=DataFrame(info), geno=SimpleList(geno))
}

.rowDataNames <- function(vcf, rowData, ...)
{
    idx <- vcf$ID == "."
    vcf$ID[idx] <- paste("chr", seqnames(rowData[idx]), ":", 
        start(rowData[idx]), sep="") 
    names(rowData) <- vcf$ID
    rowData
}


.toDNAStringSet <- function(x)
{
    xx <- unlist(strsplit(x, ",", fixed=TRUE))
    ulist <- gsub("\\.", "", xx)
    ulist[is.na(ulist)] <- ""
    DNAStringSet(ulist)
}

.toDNAStringSetList <- function(x)
{
    dna <- .toDNAStringSet(x)
    idx <- elementLengths(strsplit(x, ","))
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
    ## FIXME: should C return matrices only
    idx <- which(lapply(x, is.array) == TRUE)
    if (length(idx) != 0)
        for (i in idx) 
            x[[i]] <- DataFrame(I(matrix(x[[i]], ncol=dim(x[[i]])[3])))
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

## for PolyPhen and SIFT :

.getKcol <- function(conn)
{
    sql <- "SELECT value FROM metadata WHERE name='Key column'"
    as.character(dbGetQuery(conn, sql))
}

## convert character vector into an SQL IN condition
.sqlIn <- function(vals)
{
    if (length(vals) == 0L)
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

