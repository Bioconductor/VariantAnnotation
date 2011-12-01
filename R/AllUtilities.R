## helpers for readVcf :

.VcfToSummarizedExperiment <- function(vcf, file, genome, ...)
{
    vcf <- vcf[[1]]
    hdr <- scanVcfHeader(file)[[1]][["Header"]]

    ## assays
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

    ## rowdata
    ref <- .toDNAStringSet(vcf$REF)
    alt <- CharacterList(as.list(vcf$ALT))
    info <- .parseINFO(vcf$INFO)
    meta <- DataFrame(REF=ref, ALT=alt, QUAL=vcf$QUAL, 
                      FILTER=vcf$FILTER, info)

    rowData <- GRanges(seqnames=Rle(vcf$CHROM), 
                       ranges=IRanges(start=vcf$POS, width=width(ref)))
    values(rowData) <- meta 
    names(rowData) <- vcf$ID
    genome(seqinfo(rowData)) <- genome

    ## colData
    if (length(vcf$GENO) > 0) {
        sampleID <- colnames(vcf$GENO[[1]]) 
        samples <- length(sampleID) 
        colData <- DataFrame(
                     Samples=seq_len(samples),
                     row.names=sampleID)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    SummarizedExperiment(
      assays=geno, exptData=SimpleList(HEADER=hdr),
      colData=colData, rowData=rowData)
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
    ismissing <- grep(".", x, fixed=TRUE)
    eltlen <- sapply(as.list(x), 
        function(i) {
            length(unlist(strsplit(i, ",", fixed=TRUE)))
        })
    eltlen[ismissing] <- 0
    pbw <- PartitioningByWidth(eltlen)
    unlisted <- .toDNAStringSet(x)
    IRanges:::newCompressedList("DNAStringSetList", 
                                unlistData=unlisted[width(unlisted) > 0], 
                                end=end(pbw))
}

.parseINFO <- function(x)
{
    if (is.list(x)) {
        cmb <- lapply(x, 
            function(elt) {
                if (is.list(elt)) {
                    dat <- CharacterList(elt)
                } else { 
                    dat <- data.frame(elt)
                    names(dat) <- seq_len(ncol(dat))
                }
                dat 
            })
        info <- DataFrame(cmb) 
    } else {
        info <- DataFrame(x) 
    }
    if (is.null(names(x)))
        colnames(info) <- "INFO"
    info 
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


## helpers for PolyPhen and SIFT

.getKcol <- function(conn)
{
    sql <- "SELECT value FROM metadata WHERE name='Key column'"
    as.character(dbGetQuery(conn, sql))
}

### convert character vector into an SQL IN condition
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

duplicateRSID <- function(db, keys, ...)
{
    fmtrsid <- .sqlIn(keys)
    sql <- paste("SELECT * FROM duplicates WHERE RSID IN (",
                 fmtrsid, ")", sep="")
    q1 <- dbGetQuery(db$conn, sql)

    fmtgp <- .sqlIn(unique(q1$DUPLICATEGROUP))
    gpsql <- paste("SELECT * FROM duplicates WHERE DUPLICATEGROUP IN (",
                   fmtgp, ")", sep="")
    q2 <- dbGetQuery(db$conn, gpsql)

    matched <- q2[!q2$RSID %in% keys, ]
    matchedlst <- split(matched$RSID, matched$DUPLICATEGROUP)
    names(matchedlst) <- q1$RSID[match(names(matchedlst), q1$DUPLICATEGROUP)]

    missing <- !keys %in% q2$RSID
    if (any(missing)) {
        warning(paste("keys not found in database : ", keys[missing],
                      sep=""))
        missinglst <- list(rep(NA, sum(missing)))
        names(missinglst) <- keys[missing]
        matchedlst <- c(matchedlst, missinglst)
    }

    matchedlst[order(match(names(matchedlst), keys))]
}

