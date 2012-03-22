library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_predictCoding_empty <- function()
{
    data <- GRanges(
        seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(263892,264258, 5538757), width=1))
    alt=DNAStringSet(c("G", "T", "A"))
    res <- predictCoding(data, txdb, Hsapiens, alt)
    checkIdentical(GRanges(), res)
}

test_predictCoding_varAllele <- function()
{
    data <- GRanges(
        seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(51153371, 51153463, 51153466), width=1))
    alt=DNAStringSet(c("C", "T", ""))
    res <- suppressWarnings(predictCoding(data, txdb, Hsapiens, alt))
    expected <- DNAStringSet(c("C", "C", "T", "T", "", ""))
    checkIdentical(length(values(res)[["varAllele"]]), length(expected)) 
    checkIdentical(unlist(values(res)[["varAllele"]]), unlist(expected)) 
}

test_refLocsToLocalLocs <- function()
{
    data <- GRanges(
        seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(51153371, 51153463, 51153466), width=1))
    res <- refLocsToLocalLocs(data, txdb)
    expected <- IntegerList(737, 184, 767, 214, 768, 215)
    checkIdentical(values(res)[["protein_loc"]], expected)

    expected <- IRanges(c(2209, 550, 2301, 642, 2304, 645),
                        c(2209, 550, 2301, 642, 2304, 645))
    checkIdentical(values(res)[["cds_loc"]], expected)

    expected <- IRanges(c(2209, 832, 2301, 924, 2304, 927),
                        c(2209, 832, 2301, 924, 2304, 927))
    checkIdentical(values(res)[["cDNA_loc"]], expected)
} 

