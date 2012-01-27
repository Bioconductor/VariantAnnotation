library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

test_predictCoding_empty <- function()
{
    data <- GRanges(seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(263892,264258, 5538757), width=1))
    alt=DNAStringSet(c("G", "T", "A"))
    aa <- predictCoding(data, txdb, seqSource=Hsapiens,
        varAllele=alt)
    checkTrue(all(colnames(aa) %in%  c("queryID", "txID", "refSeq", "varSeq",
        "refAA", "varAA", "consequence", "geneID", "cdsID")))
    emptyDF <- DataFrame(queryID=integer(), consequence=character(),
        refSeq=DNAStringSet(), varSeq=DNAStringSet(),
        refAA=AAStringSet(), varAA=AAStringSet(),
        txID=integer(), geneID=character(), cdsID=integer()) 
    checkIdentical(aa, emptyDF)
}

test_predictCoding_varAllele <- function()
{
    data <- GRanges(seqnames=rep("chr22", 3),
        ranges=IRanges(start=c(51153371, 51153463, 51153466), width=1))
    alt=DNAStringSet(c("G", "T", ""))
    aa <- suppressWarnings(predictCoding(data, txdb, seqSource=Hsapiens,
        varAllele=alt))
    checkTrue(!which(width(alt) == 0) %in%
        unique(aa$queryID))
}

