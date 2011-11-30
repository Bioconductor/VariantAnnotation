library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

test_locateVariants <- function()
{
    data <- GRanges(
        seqnames=rep("chr16", 3),
        ranges=IRanges(start=c(97430, 363400, 90124759), 
                       end=c(97434, 363400, 90124760)))
    loc <- locateVariants(data, txdb)

    checkEquals(ncol(loc), 4)
    checkIdentical(colnames(loc), c("queryHits", "txID", "geneID", "Location"))
    checkTrue(all(unique(loc$queryHits) %in% seq_len(length(data))))
    checkTrue(all(unique(loc$Location) %in% c("intron", "coding", "3UTR",
        "5UTR", "transcript_region", "intergenic", "NA")))
}

