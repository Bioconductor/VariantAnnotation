library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
cdsByTx <- cdsBy(txdb, "tx") 

test_predictCoding_empty <- function()
{
    data <- 
        GRanges(
          seqnames=rep("chr22", 3),
          ranges=IRanges(start=c(263892,264258, 5538757), width=1),
          alt=DNAStringSet(c("G", "T", "A")))
    aa <- predictCoding(data, cdsByTx, seqSource=Hsapiens,
        varAllele="alt")
    checkTrue(all(colnames(aa) %in%  c("queryHits", "txID", "refSeq", "varSeq",
        "refAA", "varAA", "Consequence")))

    emptyDF <- DataFrame(queryHits=character(0), txID=character(0),
               refSeq=DNAStringSet(), varSeq=DNAStringSet(),
               refAA=AAStringSet(), varAA=AAStringSet(),
               Consequence=character(0))
    checkIdentical(aa, emptyDF)
}

test_predictCoding_varAllele <- function()
{
    data <- 
        GRanges(
          seqnames=rep("chr22", 3),
          ranges=IRanges(start=c(51153371, 51153463, 51153466), width=1),
          alt=DNAStringSet(c("G", "T", "A")))
    checkException(predictCoding(data, cdsByTx, seqSource=Hsapiens,
        varAllele="var"), silent=TRUE)

    values(data)["alt"] <- c("G", "T", "A")
    checkException(predictCoding(data, cdsByTx, seqSource=Hsapiens,
        varAllele="alt"), silent=TRUE)

    values(data)["alt"] <- DNAStringSet(c("G", "T", ""))
    aa <- predictCoding(data, cdsByTx, seqSource=Hsapiens,
        varAllele="alt")
    checkTrue(!which(width(values(data)[["alt"]]) == 0) %in%
        unique(aa$queryHits))
}


