.VcfToSummarizedExperiment <- function(vcf, header, vcfRanges=FALSE, ...)
{
    vcf <- vcf[[1]]
    # assays
    if (length(vcf$GENO) > 0)
        geno <- vcf$GENO 
    else
        geno <- list() 

    # rowdata
    lst <- as.list(vcf$ALT)
    lstSplit <- sapply(lst, function(x) length(unlist(strsplit(x, ",", fixed=TRUE))))
    lstSplit[lstSplit == 0] <- 1
    df <- DataFrame(POS=vcf$POS, QUAL=vcf$QUAL, FILTER=vcf$FILTER)
    info <- as.data.frame(vcf$INFO) 
    df <- append(df, info)

    reference <- gsub("\\.", "", vcf$REF)
    reference[is.na(reference)] <- ""
    reference <- DNAStringSet(reference)

    alt <- unlist(strsplit(vcf$ALT, ",", fixed=TRUE))
    alt <- gsub("\\.", "", alt)
    alt[is.na(alt)] <- ""
    alt <- DNAStringSet(alt)
    if (vcfRanges)  {
        allalt <- DataFrame(alt)
        altDF <- split(allalt, rep(seq_len(length(lstSplit)), lstSplit))
    } else { 
        altDF <- .convertToBiocRanges(reference, alt, vcf, lstSplit)
    }

    rowdata <- GRanges(Rle(vcf$CHROM), IRanges(vcf$POS, width=1))
    df$REF <- reference
    df$ALT <- altDF 
    values(rowdata) <- df
    names(rowdata) <- vcf$ID

    ## colData
    if (length(vcf$GENO) > 0) {
        sampleID <- colnames(vcf$GENO[[1]]) 
        samples <- length(sampleID) 
        coldata <- DataFrame(Samples=seq_len(samples),
            row.names=sampleID)
    } else {
    coldata <- DataFrame(Samples=character(0))
    }

    SummarizedExperiment(assays=geno, colData=coldata, rowData=rowdata)
}

.convertToBiocRanges <- function(reference, alt, vcf, lstSplit, ...)
{
    ## snp
    ref <- rep(reference, lstSplit)
    start <- end <- rep(vcf$POS, lstSplit)
    ## deletion
    deletion <- width(ref) > width(alt)
    start[deletion] <- start[deletion] + 1
    end[deletion] <- start[deletion] +
        abs(width(ref) - width(alt))[deletion]
    ## insertion
    insertion <- width(ref) < width(alt)
    start[insertion] <- start[insertion] + 2
    end[insertion] <- end[insertion] + 1
    ## block substitution 
    blockSub <- (width(ref) == width(alt)) & (width(ref) != 1)
    start[blockSub] <- start[blockSub] + 1
    end[blockSub] <- start[blockSub] + width(ref)[blockSub] - 1

    ## remove leading allele in all except snps and monomorphic 
    notSNPorMONO_ext <- width(ref) > 1
    notSNPorMONO_con <- width(reference) > 1
    reference[notSNPorMONO_con] <-
        narrow(reference[notSNPorMONO_con], start=2)
    alt[notSNPorMONO_ext] <- narrow(alt[notSNPorMONO_ext], start=2)
    res <- DataFrame(alt, ranges=IRanges(start, end))
    split(res, rep(seq_len(length(lstSplit)), lstSplit))
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


