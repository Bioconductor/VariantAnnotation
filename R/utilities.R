.VcfToSummarizedExperiment <- function(vcf, file, ...)
{
    vcf <- vcf[[1]]

    ## assays
    if (length(vcf$GENO) > 0) {
       geno <- lapply(vcf$GENO, function(elt) {
           if (is.list(elt))
               matrix(unlist(elt, recursive=FALSE), ncol=1) 
           else
               elt
       })
    } else {
        geno <- list() 
    }

    ## rowdata
    lst <- as.list(vcf$ALT)
    ismissing <- grep(".", lst, fixed=TRUE)
    eltlen <- 
        sapply(lst, function(x) {
            length(unlist(strsplit(x, ",", fixed=TRUE)))
        })
    #eltlen[eltlen == 0] <- 1
    eltsplit <- rep(seq_len(length(eltlen)), eltlen)

    ref <- .toDNAStringSet(vcf$REF)
    alt <- DataFrame(ALT=.toDNAStringSet(vcf$ALT))
    altsplit <-  split(alt, eltsplit)
    if (any(ismissing))
        altsplit[[ismissing]] <- DataFrame() 
    info <- data.frame(vcf$INFO) 
    if (is.null(names(vcf$INFO)))
        colnames(info) <- "INFO"
    DF <- DataFrame(REF=ref, ALT=NA, QUAL=vcf$QUAL, 
                    FILTER=vcf$FILTER, data.frame(info))
    DF$ALT <- altsplit 

    rowData <- GRanges(Rle(vcf$CHROM), 
                       IRanges(start=vcf$POS, 
                       width=width(ref)))
    values(rowData) <- DF 
    names(rowData) <- vcf$ID

    ## colData
    if (length(vcf$GENO) > 0) {
        sampleID <- colnames(vcf$GENO[[1]]) 
        samples <- length(sampleID) 
        colData <- DataFrame(Samples=seq_len(samples),
            row.names=sampleID)
    } else {
        colData <- DataFrame(Samples=character(0))
    }

    ## exptData
    header <- scanVcfHeader(file)[[1]][["Header"]]

    SummarizedExperiment(assays=geno, exptData=SimpleList(HEADER=header),
                         colData=colData, rowData=rowData)
}

.toDNAStringSet <- function(x)
{
    ulist <- unlist(strsplit(x, ",", fixed=TRUE))
    ulist <- gsub("\\.", "", ulist)
    ulist[is.na(ulist)] <- ""
    DNAStringSet(ulist)
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

