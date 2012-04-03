library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 

cdsbytx <- cdsBy(txdb)
intbytx <- intronsByTranscript(txdb)
txbygene <- transcriptsBy(txdb, "gene")
gr <- GRanges("chr22", 
              IRanges(c(16268137, 16287254, 16190792, 16164570,
                        18209442, 18121652, 24314750, 25508661), 
                      width=c(1,1,1,1,3,3,2,2)),
              strand=c("-", "-", "-", "+", "+", "+", "+", "+"))
 
test_locateVariants_subject <- function()
{
    cols <- c("location", "queryID", "txID", "cdsID")
    loc1 <- locateVariants(gr, txdb, CodingVariants())
    loc2 <- locateVariants(gr, cdsbytx, CodingVariants())
    checkIdentical(values(loc1)[ ,cols], values(loc2)[ ,cols])
 
    loc1 <- locateVariants(gr, txdb, IntronVariants())
    loc2 <- locateVariants(gr, intbytx, IntronVariants())
    checkIdentical(values(loc1)[ ,cols], values(loc2)[ ,cols])

    loc1 <- locateVariants(gr, txdb, SpliceSiteVariants())
    loc2 <- locateVariants(gr, intbytx, SpliceSiteVariants())
    checkIdentical(values(loc1)[ ,cols], values(loc2)[ ,cols])
 
    loc1 <- locateVariants(gr, txdb, IntergenicVariants())
    loc2 <- locateVariants(gr, txbygene, IntergenicVariants())
    checkIdentical(values(loc1)[ ,cols], values(loc2)[ ,cols])
}

test_locateVariants_query <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    vcf <- renameSeqlevels(vcf, c("20" = "chr20"))
    loc1 <- locateVariants(vcf, txdb, IntergenicVariants())
    loc2 <- locateVariants(rowData(vcf), txdb, IntergenicVariants())
    checkIdentical(loc1, loc2) 
}

