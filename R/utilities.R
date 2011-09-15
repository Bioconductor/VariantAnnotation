.VcfToSummarizedExperiment <- function(vcf, ...)
{
    vcf <- vcf[[1]]
    ## assays
    if (length(vcf$GENO) > 0)
        geno <- vcf$GENO 
    else
        geno <- list() 

    ## rowdata
    lst <- as.list(vcf$ALT)
    lstSplit <- sapply(lst, function(x) length(unlist(strsplit(x, ",", fixed=TRUE))))
    lstSplit[lstSplit == 0] <- 1
    ref <- .toDNAStringSet(vcf$REF)
    alt <- DataFrame(.toDNAStringSet(vcf$ALT))
    df <- DataFrame(REF=ref, ALT=NA, QUAL=vcf$QUAL, 
        FILTER=vcf$FILTER, data.frame(vcf$INFO))
    df$ALT <- split(alt, rep(seq_len(length(lstSplit)), lstSplit))

    rowData <- GRanges(Rle(vcf$CHROM), 
        IRanges(start=vcf$POS, width=width(ref)))
    values(rowData) <- df
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

    SummarizedExperiment(assays=geno, colData=colData, rowData=rowData)
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


