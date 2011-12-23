### =========================================================================
### Helper functions not exported 
### =========================================================================

## for readVcf :

.VcfToSummarizedExperiment <- function(vcf, file, genome, ..., param)
{
    vcf <- vcf[[1]]
    hdr <- scanVcfHeader(file)[[1]][["Header"]]

    ## assays
    if (length(vcf$GENO) > 0) {
        gno <- lapply(vcf$GENO, 
            function(elt) {
                if (is.list(elt))
                    do.call(rbind, elt) 
                else
                    elt
            })
    } else {
        gno <- list() 
    }

    ## info 
    if (!is.null(param)) {
        if (class(param) == "ScanVcfParam") {
            if (!identical(character(), vcfInfo(param)))
                vcf$INFO <- vcf$INFO[names(vcf$INFO) %in% vcfInfo(param)]
        }
    } 
    inf <- .parseINFO(vcf$INFO,  hdr[["INFO"]])

    ## rowdata
    ref <- .toDNAStringSet(vcf$REF)
    alt <- CharacterList(as.list(vcf$ALT))
    meta <- DataFrame(REF=ref, ALT=alt, QUAL=vcf$QUAL, FILTER=vcf$FILTER)
    rowData <- GRanges(seqnames=Rle(vcf$CHROM), 
                       ranges=IRanges(start=vcf$POS, width=width(ref)))
    values(rowData) <- meta
    idx <- vcf$ID == "."
    vcf$ID[idx] <- paste("chr", seqnames(rowData[idx]), ":", start(rowData[idx]), sep="") 
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

    VCF(assays=SimpleList(gno), colData=colData, rowData=rowData, 
        exptData=SimpleList(HEADER=hdr), info=SimpleList(inf))
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
           x[[i]] <- 
             switch(type[[i]],
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

