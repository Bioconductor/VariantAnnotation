library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_predictCoding_empty <- function()
{
    data <- GRanges(
        seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(263892,264258, 5538757), width=1))
    alt <- DNAStringSet(c("G", "T", "A"))
    current <- predictCoding(data, txdb, Hsapiens, alt)
    expected <- GRanges()
    checkIdentical(current, expected)
}

test_predictCoding_varAllele <- function()
{
    data <- GRanges(
        seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(51153371, 51153463, 51153466), width=1))
    alt <- DNAStringSet(c("C", "T", ""))
    current <- suppressWarnings(predictCoding(data, txdb, Hsapiens, alt))
    expected <- DNAStringSet(c("C", "C", "T", "T", "", ""))
    checkIdentical(length(values(current)[["varAllele"]]), length(expected)) 
    checkIdentical(unlist(values(current)[["varAllele"]]), unlist(expected)) 
}

test_refLocsToLocalLocs <- function()
{
    data <- GRanges(
        seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(51153371, 51153463, 51153466), width=1))
    current <- refLocsToLocalLocs(data, txdb)
    expected <- IntegerList(737, 184, 767, 214, 768, 215)
    checkIdentical(values(current)[["proteinLoc"]], expected)

    expected <- IRanges(c(2209, 550, 2301, 642, 2304, 645),
                        c(2209, 550, 2301, 642, 2304, 645))
    checkIdentical(values(current)[["cdsLoc"]], expected)
} 

