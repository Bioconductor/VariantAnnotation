.parseTabix <- function(tbx, header, param, ...)
{
    ## currently return all fields in ScanBcfParam
    tmpl <- Rsamtools:::.bcf_template(param=ScanBcfParam())
    tmpl$RecordsPerRange <- sum(elementLengths(tbx)) 
    strings <- strsplit(unlist(tbx, use.names=FALSE), "\t", fixed=TRUE)
    nfields <- length(strings[[1]]) 
    mat <- matrix(unlist(strings, use.names=FALSE), ncol=nfields, byrow=TRUE)
    if (!is.null(header[[1]]$Sample))
        nsamples <- length(header[[1]]$Sample)
    else 
        nsamples <- 0

    ## required fields 
    tmpl[1:8] <- lapply(seq_len(8), function(x, mat) mat[,x], mat)
    tmpl[["POS"]] <- as.numeric(tmpl[["POS"]])
    tmpl[["QUAL"]] <- gsub(".", "", tmpl[["QUAL"]])
    tmpl[["QUAL"]] <- as.numeric(tmpl[["QUAL"]])
    #if (all(tmpl[["INFO"]] == ".")) 
    #    tmpl[["INFO"]] <- rep(NA, length(tmpl[["INFO"]]))

    ## optional fields 
    if (nsamples > 0) {
        tmpl[["FORMAT"]] <- mat[,9]
        samples <- mat[,-c(1:9)]
        formats <- unlist(strsplit(tmpl[["FORMAT"]][1], ":", fixed=TRUE))
        tmpl[["FORMAT"]] <- formats
        tmpl[["GENO"]][formats] <- .parseGENO(tmpl, formats, samples, header)
        tmpl[["GENO"]] <-  tmpl[["GENO"]][formats]
    } else {
        tmpl[["GENO"]] <- list()
    }

    tmpl
}


.VcfToSummarizedExperiment <- function(vcf, header, raw=FALSE, ...)
{
    # assays
    if (length(vcf$GENO) > 0)
        geno <- lapply(vcf$GENO, function(x) matrix(x, 
           nrow=vcf$RecordsPerRange, byrow=FALSE))
    else
        geno <- list() 

    # rowdata
    lst <- as.list(vcf$ALT)
    lstSplit <- sapply(lst, function(x) length(unlist(strsplit(x, ",", fixed=TRUE))))
    lstSplit[lstSplit == 0] <- 1
    df <- DataFrame(POS=vcf$POS, QUAL=vcf$QUAL, FILTER=vcf$FILTER)
    info <- .parseINFO(vcf$INFO, header)
    df <- append(df, info)

    reference <- gsub("\\.", "", vcf$REF)
    reference[is.na(reference)] <- ""
    reference <- DNAStringSet(reference)

    alt <- unlist(strsplit(vcf$ALT, ",", fixed=TRUE))
    alt <- gsub("\\.", "", alt)
    alt[is.na(alt)] <- ""
    alt <- DNAStringSet(alt)
    if (raw)  {
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
        sampleID <- header[[1]]$Sample
        samples <- length(sampleID) 
        coldata <- DataFrame(Samples=seq_len(samples),
            row.names=sampleID)
    } else {
    coldata <- DataFrame(Samples=character(0))
    }

    SummarizedExperiment(assays=geno, colData=coldata, rowData=rowdata)
}


.vcfDataTypeConversion <- function(x, ...)
{
    x[x == "String"] <- "character"
    x[x == "Integer"] <- "integer"
    x[x == "Flag"] <- "factor"
    x[x == "Float"] <- "numeric"
    x
}

.parseINFO <- function(rawinfo, header, ...)
{
    ## FIXME : handle case of > 1 value numeric, integer 
    ## FIXME : force type=flag to factor 
    if (all(rawinfo == ".")) 
        return(DataFrame())

    info <- strsplit(rawinfo, ";", fixed=TRUE)
    u <-  do.call(rbind, strsplit(unlist(info), "=", fixed=TRUE))
    u <- u[u[,1] != ".", ]
    keys <- u[,1]
    values <- u[,2]
    infoHeader <- header[[1]]$Header$INFO
    flags <- rownames(infoHeader)[infoHeader$Type == "Flag"]
    values[values %in% flags] <- 1
    lstvalues <- split(values, keys)
    uniquekeys <- names(lstvalues) 
    lstkeys <-  lapply(uniquekeys, function(x, info) grep(x, info), info)

    vcfType <- infoHeader$Type[match(uniquekeys, rownames(infoHeader))]
    RType <- .vcfDataTypeConversion(vcfType)

    lst <- lapply(seq_len(length(lstkeys)), 
        function(x, lstkeys, lstvalues, RType) {
            new <- numeric(length(info))
            new[lstkeys[[x]]] <- switch(RType[x], 
                character = as.character(lstvalues[[x]]),
                numeric = as.numeric(lstvalues[[x]]),
                integer = as.integer(lstvalues[[x]]),
                factor = as.factor(lstvalues[[x]]))
            new
        }, lstkeys, lstvalues, RType)

    DF <- as(lst, "DataFrame")
    colnames(DF) <- uniquekeys
    DF
}

.parseGENO <- function(tmpl, formats, rawgeno, header, ...)
{
    ## FIXME : handle case of > 1 value for numeric, integer 
    ## FIXME : force type=flag to factor 
    strings <- strsplit(unlist(rawgeno, use.names=FALSE), ":", fixed=TRUE)
    formatHeader <- header[[1]]$Header$FORMAT[formats %in%
        rownames(header[[1]]$Header$FORMAT)]
    vcfType <- formatHeader$Type[match(formats, rownames(formatHeader))]
    formatType <- .vcfDataTypeConversion(vcfType)
    mat <- matrix(unlist(strings, use.names=FALSE), ncol=length(formats), byrow=TRUE)
    lst <- lapply(seq_len(length(formatType)), function(x, mat, formatType){
            new <- mat[,x]
            #mode(new) <- formatType[x]
            matrix(new, ncol=ncol(rawgeno))
        }, mat, formatType)
    names(lst) <- formats
    lst
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


