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
        ranges=IRanges(start=c(51153371, 51153462, 51153466), ## altpos = 1,2,3
        width=rep(1, 3)))

    ## missing alt allele
    alt <- DNAStringSet(c("G", "T", ""))
    current <- suppressWarnings(predictCoding(data, txdb, Hsapiens, alt))
    expected <- DNAStringSet(c("G", "G", "T", "T", "", ""))
    checkIdentical(unlist(values(current)[["varAllele"]]), unlist(expected))

    ## width = 1
    checkEquals(width(values(current)[["refSeq"]]), rep(3L, 6))

    ## width = 2
    width(data) <- rep(2, 3) 
    alt <- DNAStringSet(c("GG", "TT", ""))
    current <- suppressWarnings(predictCoding(data, txdb, Hsapiens, alt))
    checkEquals(width(values(current)[["refSeq"]]), c(rep(3L, 4), 6L, 6L))
 
    ## width = 3 
    width(data) <- rep(3, 3) 
    alt <- DNAStringSet(c("GGG", "TTT", ""))
    current <- suppressWarnings(predictCoding(data, txdb, Hsapiens, alt))
    checkEquals(width(values(current)[["refSeq"]]), c(3L, 3L, rep(6L, 4))) 
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

