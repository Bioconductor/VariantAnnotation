library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(BSgenome.Dmelanogaster.UCSC.dm3)
txdb <- Dmelanogaster_UCSC_dm3_ensGene_TxDb


test_predictCoding <- function()
{
    data <- GRanges(seqnames=c("chr4", "chr4", "chr2L"),
                IRanges(start=c(263892,264258, 5538757), width=1),
                alt=DNAStringSet(c("G", "T", "A")))
    aa <- predictCoding(data, txdb, seqSource=Dmelanogaster,
        varAllele="alt")

    checkEquals(ncol(aa), 10)
    checkIdentical(colnames(aa), c("queryHits", "txID", "refSeq", "varSeq",
        "refAA", "varAA", "Consequence", "cds_id", "cds_name", "exon_rank"))
}

test_predictCoding_no_match <- function()
{
    data <- GRanges(seqnames=c("chr4", "chr2L"),
                IRanges(start=c(264258, 5538750), width=1),
                alt=DNAStringSet(c("T", "A")))
    aa <- predictCoding(data, txdb, seqSource=Dmelanogaster,
        varAllele="alt")
    emptyDF <- DataFrame(queryHits=character(0), txID=character(0),
               refSeq=DNAStringSet(), varSeq=DNAStringSet(),
               refAA=AAStringSet(), varAA=AAStringSet(),
               Consequence=character(0))
 
    checkIdentical(aa, emptyDF)
}

test_predictCoding_varAllele <- function()
{
    data <- GRanges(seqnames=c("chr4", "chr2L"),
                IRanges(start=c(264258, 5538750), width=1),
                alt=DNAStringSet(c("T", "A")))

    checkException(predictCoding(data, txdb, seqSource=Dmelanogaster,
        varAllele="var"), silent=TRUE)
    checkException(predictCoding(data, txdb, seqSource=Dmelanogaster))
}

