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
    cols <- c("LOCATION", "QUERYID", "TXID", "CDSID")
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

test_locateVariants_strand <- function()
{
    gr <- GRanges("chr1", IRanges(c(12190, 12595, 13403), width=1), "-")
    loc1 <- locateVariants(gr, cdsbytx, CodingVariants(), 
                           ignore.strand=TRUE)
    checkIdentical(c(1L, 2L, 3L), values(loc1)[["QUERYID"]]) 
    loc2 <- locateVariants(gr, cdsbytx, CodingVariants(), 
                           ignore.strand=FALSE)
    checkIdentical(integer(), values(loc2)[["QUERYID"]]) 
    loc1 <- locateVariants(gr, cdsbytx, SpliceSiteVariants(), 
                           ignore.strand=TRUE)
    checkIdentical(c(1L, 2L, 3L), values(loc1)[["QUERYID"]]) 
    loc2 <- locateVariants(gr, cdsbytx, SpliceSiteVariants(), 
                           ignore.strand=FALSE)
    checkIdentical(integer(), values(loc2)[["QUERYID"]]) 
}

.extract <- function(x, col) as.vector(values(x)[[col]])
test_FlankingVariants <- function()
{
    ## empty 
    q <- GRanges("chr1", IRanges(15, width=1), "+")
    s <- GRanges(c("chr1", "chr1"), IRanges(c(10, 30), width=11), "+")
    current <- locateVariants(q, s, FlankingVariants(5, 5))
    checkTrue(length(current) == 0)
 
    ## endpoint
    q <- GRanges("chr1", IRanges(20, width=1), "+")
    s <- GRanges(c("chr1", "chr1"), IRanges(c(10, 30), width=11), "+")
    current <- locateVariants(q, s, FlankingVariants(5, 5))
    checkTrue(length(current) == 0)
 
    q <- GRanges(c("chr1", "chr1"), IRanges(c(8, 42), width=1), "+")
    s <- GRanges(c("chr1", "chr1"), IRanges(c(10, 30), width=11), "+")
    current <- locateVariants(q, s, FlankingVariants(5, 5))
    target <- c("upstream", "downstream")
    checkEquals(target, .extract(current, "LOCATION"))

    q <- GRanges(rep("chr2", 4), IRanges(c(10, 8, 42, 50), width=1), "+")
    s <- GRanges(c("chr2", "chr2"), IRanges(c(10, 30), width=11), "+")
    current <- locateVariants(q, s, FlankingVariants(5, 5))
    target <- c("upstream", "downstream")
    checkEquals(target, .extract(current, "LOCATION"))
    target <- c(2, 3)
    checkEquals(target, .extract(current, "QUERYID"))
}
