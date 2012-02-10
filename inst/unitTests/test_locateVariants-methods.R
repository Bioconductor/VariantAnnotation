library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

test_locateVariants <- function()
{
    data <- GRanges(seqnames=rep("chr16", 3), 
        ranges=IRanges(start=c(97430, 363400, 90124759), 
        end=c(97434, 363400, 90124760)))
    loc <- locateVariants(data, txdb)

    checkIdentical(colnames(loc), 
        c("queryID", "location", "txID", "cdsID", "geneID", "precedesID",
          "followsID"))
    checkTrue(all(unique(loc$queryID) %in% seq_len(length(data))))
    checkTrue(all(unique(loc$location) %in% c("intron", "coding", "3UTR",
        "5UTR", "transcript_region", "intergenic", "NA")))
}

