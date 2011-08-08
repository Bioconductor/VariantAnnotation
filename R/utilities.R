## TODO : - Modify scanTabix/scanBcf to pull over 
##          HQ and FT from INFO and sample names
##        - Make scanTabix work with ScanBcfParam 

.parseTabix <- function(tbx, param, ...)
{
    ## currently returning all fields in ScanBcfParam
    tmpl <- Rsamtools:::.bcf_template(param=ScanBcfParam())
    ulst <- unlist(tbx, use.names=FALSE)
    strings <- strsplit(ulst, "\t", fixed=TRUE)
    elt <- elementLengths(tbx)
    nfields <- length(strings[[1]]) 
    tmpl$RecordsPerRange <- sum(elt) 

    mat <- matrix(unlist(strings, use.names=FALSE), ncol=nfields, byrow=TRUE)

    ## 8 required fields
    reqfields <- 8
    if (nfields > reqfields) {
        nsamples <- nfields - (reqfields + 1)
        fill <- reqfields + 1
    } else {
        nsamples <- 0
        fill <- reqfields
    }

    ## FIXME : datatypes
    tmpl[1:fill] <- lapply(seq_len(fill), function(x, mat) mat[,x], mat)
    tmpl[["POS"]] <- as.numeric(tmpl[["POS"]])

    ## geno
    if (nfields > reqfields) {
        nsamples <- nfields - (reqfields + 1) 
        samples <- mat[,c((reqfields + 2):nfields)]
        idx <- unlist(strsplit(tmpl[["FORMAT"]][1], ":"))
        ss <- strsplit(samples, ":", fixed=TRUE)
        genomat <- matrix(unlist(ss, use.names=FALSE), ncol=length(idx), byrow=TRUE)
        genolst <- lapply(seq_len(length(idx)), function(x, genomat){
            matrix(genomat[,x], ncol=nsamples)
        }, genomat)
        names(genolst) <- idx
        tmpl[["GENO"]][idx] <- genolst[idx]
        tmpl[["GENO"]] <- tmpl[["GENO"]][idx]
    } else { 
        tmpl <- tmpl[1:fill] 
    }

    tmpl
}


.VcfToSummarizedExperiment <- function(vcf, raw=FALSE, ...)
{
    # assays
    if (!is.null(vcf$GENO))
        geno <- lapply(vcf$GENO, as.matrix)
    else
        geno <- list() 

    # rowdata
    lst <- as.list(vcf$ALT)
    lstSplit <- sapply(lst, function(x) length(unlist(strsplit(x, ","))))
    lstSplit[lstSplit == 0] <- 1
    df <- DataFrame(pos=vcf$POS, qual=vcf$QUAL, filter=vcf$FILTER,
        info=vcf$INFO)
    reference <- DNAStringSet(gsub("\\.", "", vcf$REF))
    alt <- unlist(strsplit(vcf$ALT, ","))
    alt <- DNAStringSet(gsub("\\.", "", alt))

    if (raw) {
        alt <- DataFrame(alt)
    } else {
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
        alt <- DataFrame(alt, ranges=IRanges(start, end))
    }

    df$ref <- reference
    df$alt <- split(alt, rep(seq_len(length(lstSplit)), lstSplit))
    rowdata <- GRanges(seqnames=Rle(vcf$CHROM),
                       ranges=IRanges(start=vcf$POS, width=1))
    values(rowdata) <- df
    names(rowdata) <- vcf$ID

    ## colData
    if (!is.null(vcf$GENO)) {
        samples <- ncol(do.call(rbind, geno))
        coldata <- DataFrame(Samples=seq_len(samples),
            row.names=LETTERS[1:samples])
       } else coldata <- DataFrame(Samples=1, row.names=LETTERS[1])

    SummarizedExperiment(assays=geno, colData=coldata, rowData=rowdata)
}


.cumsumShifted <- function(x) 
{
    cs <- unlist(cumsum(x))
    shifted <- c(0L, head(cs, -1))
    shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
    shifted
}


