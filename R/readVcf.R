## TODO :
## imprecise variants ALT field may contain descriptive alleles
## imprecise should be marked by imprecise flag in INFO field
## strand bias in INFO$sb

readVcf <- function(file, index=file, param=ScanBcfParam, raw=FALSE, ...)
{
    ## FIXME : change index when tabix is in place
    vcf <- scanBcf(file, index=character(0))

    # assays 
    if (length(vcf$GENO) > 0) 
        geno <- lapply(vcf$GENO, as.matrix) 
 
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
        reference[notSNPorMONO_con] <- narrow(reference[notSNPorMONO_con], start=2)
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
    ## FIXME : modify scanBcf to pull over this infomation 
    if (length(vcf$GENO) > 0) {
        samples <- ncol(do.call(rbind, geno))
        coldata <- DataFrame(Samples=seq_len(samples),
            row.names=LETTERS[1:samples])
       } else coldata <- DataFrame(Samples=1, row.names=LETTERS[1])

    SummarizedExperiment(assays=geno, colData=coldata, rowData=rowdata)
}

