### =========================================================================
### Helper functions not exported 
### =========================================================================

## for readVcf :

.VcfAsGRanges <- function(vcf, file, genome, param, ...)
{
    vcf <- vcf[[1]]
    if (vcfAsGRanges(param) == "info")
        dat <- vcf$INFO
    else if (vcfAsGRanges(param) == "geno")
        dat <- vcf$GENO
    else
        stop("asGRanges must be specified as 'info' or 'geno")
 
    rowData <- .rowData(vcf, genome)
    nvar <- length(rowData)
    idx <- lapply(dat, function(elt, nvar) {
               if (is.list(elt))
                   elementLengths(elt)
               else
                   rep(1, nvar)
           }, nvar) 
    reps <- apply(do.call(cbind, idx), 1, max)

    lst <- lapply(dat, function(elt, reps) {
               if (is(elt, "array")) {
                   slen <- rep(seq_len(length(reps)), reps)
                   data.frame(elt)[slen, ]
               } else if (is(elt, "list")) {
                   mt <- elementLengths(elt) == reps
                   newrep <- rep(1, length(mt))
                   newrep[mt == FALSE] <- reps[mt == FALSE]
                   unlist(rep(elt, newrep)) 
               } else {
                   rep(elt, reps)
               }
           }, reps)

    DF <- DataFrame(do.call(cbind, lst))
    colnames(DF) <- gsub(".X1", "", colnames(DF), fixed=TRUE)
    gr <- rowData[rep(seq_len(length(rowData)), reps)]
    values(gr) <- DF
    gr 

}

.VcfAsSummarizedExperiment <- function(vcf, file, genome, ...)
{
    vcf <- vcf[[1]]
    hdr <- scanVcfHeader(file)[[1]][["Header"]]

    ## geno 
    if (length(vcf$GENO) > 0) {
        geno <- lapply(vcf$GENO, 
            function(elt) {
                if (is.list(elt))
                    do.call(rbind, elt) 
                else
                    elt
            })
    } else {
        geno <- list() 
    }

    ## info 
    if (length(vcf$INFO) > 0)
        info <- .parseINFO(vcf$INFO,  hdr[["INFO"]])
    else 
        info <- list() 

    ## rowData
    rowData <- .rowData(vcf, genome)

    ## colData
    if (length(vcf$GENO) > 0) {
        sampleID <- colnames(geno[[1]]) 
        samples <- length(sampleID) 
        colData <- DataFrame(
                     Samples=seq_len(samples),
                     row.names=sampleID)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    VCF(geno=SimpleList(geno), colData=colData, rowData=rowData, 
        exptData=SimpleList(HEADER=hdr), info=SimpleList(info))
}

.rowData <- function(vcf, genome, ...)
{
    ref <- .toDNAStringSet(vcf$REF)
    ## FIXME: DNAStringSet when not structural
    alt <- CharacterList(as.list(vcf$ALT))
    meta <- DataFrame(REF=ref, ALT=alt, QUAL=vcf$QUAL, FILTER=vcf$FILTER)
    rowData <- GRanges(seqnames=Rle(vcf$CHROM), 
                       ranges=IRanges(start=vcf$POS, width=width(ref)))
    values(rowData) <- meta
    idx <- vcf$ID == "."
    vcf$ID[idx] <- paste("chr", seqnames(rowData[idx]), ":", start(rowData[idx]), sep="") 
    names(rowData) <- vcf$ID
    genome(seqinfo(rowData)) <- genome
    rowData
}


.toDNAStringSet <- function(x)
{
    xx <- unlist(strsplit(x, ",", fixed=TRUE))
    ulist <- gsub("\\.", "", xx)
    ulist[is.na(ulist)] <- ""
    DNAStringSet(ulist)
}

.newCompressedList <- IRanges:::newCompressedList
.parseINFO <- function(x, hdr)
{
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
    x
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

