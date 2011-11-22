library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(BSgenome.Dmelanogaster.UCSC.dm3)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene 


test_locateVariants <- function()
{
    data <- GRanges(seqnames=c("chr4", "chr4", "chr2L"),
                IRanges(start=c(263892,264258, 5538757), width=1),
                alt=DNAStringSet(c("G", "T", "A")))
    loc <- locateVariants(data, txdb)

    checkEquals(ncol(loc), 4)
    checkIdentical(colnames(loc), c("queryHits", "txID", "geneID", "Location"))
}


#test_locateVariants_unknown <- function()
#{
#    data <- GRanges(seqnames=c("chr4", "chr4"),
#                IRanges(start=c(24050, 12), width=1),
#                alt=DNAStringSet(c("T", "A")))
#    loc <- locateVariants(data, txdb)
#    DF <- DataFrame(queryHits=seq_len(length(data)),
#        txID=character(length(data)),
#        geneID=CharacterList(character(length(data))),
#        Location=factor(rep("unknown", length(data)))) 
#
#    checkIdentical(loc, DF)
#    checkEquals(levels(DF$Location), "unknown")
#}
