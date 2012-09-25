library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
cdsbytx <- cdsBy(txdb)
intbytx <- intronsByTranscript(txdb)
txbygene <- transcriptsBy(txdb, "gene")

gr <- GRanges("chr22", IRanges(c(16268137, 16287254, 16190792, 16164570,
    18209442, 18121652, 24314750, 25508661), width=c(1,1,1,1,3,3,2,2)),
    strand=c("-", "-", "-", "+", "+", "+", "+", "+"))
 
test_locateVariants_subject <- function()
{
    cols <- c("LOCATION", "QUERYID", "TXID", "CDSID")
    loc1 <- locateVariants(gr, txdb, CodingVariants())
    loc2 <- locateVariants(gr, cdsbytx, CodingVariants())
    checkIdentical(mcols(loc1)[ ,cols], mcols(loc2)[ ,cols])
 
    loc1 <- locateVariants(gr, txdb, IntronVariants())
    loc2 <- locateVariants(gr, intbytx, IntronVariants())
    checkIdentical(mcols(loc1)[ ,cols], mcols(loc2)[ ,cols])

    loc1 <- locateVariants(gr, txdb, SpliceSiteVariants())
    loc2 <- locateVariants(gr, intbytx, SpliceSiteVariants())
    checkIdentical(mcols(loc1)[ ,cols], mcols(loc2)[ ,cols])
 
    loc1 <- locateVariants(gr, txdb, IntergenicVariants())
    loc2 <- locateVariants(gr, txbygene, IntergenicVariants())
    checkIdentical(mcols(loc1)[ ,cols], mcols(loc2)[ ,cols])
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

test_locateVariants_ignore.strand <- function()
{
    gr <- GRanges("chr1", IRanges(c(12190, 12595, 13403), width=1), "-")
    loc1 <- locateVariants(gr, cdsbytx, CodingVariants(), 
                           ignore.strand=TRUE)
    checkIdentical(c(1L, 2L, 3L), mcols(loc1)$QUERYID) 
    loc2 <- locateVariants(gr, cdsbytx, CodingVariants(), 
                           ignore.strand=FALSE)
    checkIdentical(integer(), mcols(loc2)$QUERYID) 
    loc1 <- locateVariants(gr, cdsbytx, SpliceSiteVariants(), 
                           ignore.strand=TRUE)
    checkIdentical(c(1L, 2L, 3L), mcols(loc1)$QUERYID) 
    loc2 <- locateVariants(gr, cdsbytx, SpliceSiteVariants(), 
                           ignore.strand=FALSE)
    checkIdentical(integer(), mcols(loc2)$QUERYID) 
}

test_locateVariants_asHits <- function()
{
    gr <- GRanges("chr1", IRanges(c(12190, 69091, 13403), width=1))
    loc <- locateVariants(gr, cdsbytx, CodingVariants())
    hit <- locateVariants(gr, cdsbytx, CodingVariants(), asHits=TRUE)
    ## annotation element 
    loc_nms <- as.character(mcols(loc)$TXID)
    hit_nms <- names(cdsbytx[subjectHits(hit)])
    checkIdentical(loc_nms, hit_nms) 

    ## Hits lengths
    checkIdentical(length(gr), queryLength(hit))
    checkIdentical(length(cdsbytx), subjectLength(hit))
}

.extract <- function(x, col) as.vector(mcols(x)[[col]])
test_locateVariants_PromoterVariants <- function()
{
    s <- GRangesList(GRanges("chr1", IRanges(10, width=11), "+"),
                     GRanges("chr1", IRanges(30, width=11) , "+"))
    ## empty 
    q <- GRanges("chr1", IRanges(15, width=1), "+")
    current <- locateVariants(q, s, PromoterVariants(5, 5))
    checkTrue(length(current) == 0)
 
    ## endpoint
    q <- GRanges("chr1", IRanges(20, width=1), "+")
    current <- locateVariants(q, s, PromoterVariants(5, 5))
    checkTrue(length(current) == 0)
 
    ## strand 
    q <- GRanges(c("chr1", "chr1"), IRanges(c(8, 12), width=1), "+")
    current <- locateVariants(q, s, PromoterVariants(5, 5))
    checkEquals(c(1L, 2L), .extract(current, "QUERYID"))
    strand(s) <- RleList(Rle(factor("*")), Rle(factor("*")))
    strand(q) <- "*"
    current <- suppressWarnings(locateVariants(q, s, PromoterVariants(5, 5)))
    checkEquals(c(1L, 2L), .extract(current, "QUERYID"))
    q <- GRanges(c("chr1", "chr1"), IRanges(c(21, 41), width=1), "-")
    strand(s) <- RleList(Rle(factor("-")), Rle(factor("-")))
    current <- locateVariants(q, s, PromoterVariants(5, 5))
    checkEquals(c(1L, 2L), .extract(current, "QUERYID"))

    q <- GRanges(c("chr2", "chr2"), IRanges(c(9, 10), width=1), "+")
    s <- GRangesList(GRanges("chr2", IRanges(10, width=11), "+"))
    current <- locateVariants(q, s, PromoterVariants(5, 0))
    checkEquals(1L, .extract(current, "QUERYID"))
    current <- locateVariants(q, s, PromoterVariants(5, 1))
    checkEquals(c(1L, 2L), .extract(current, "QUERYID"))
    current <- locateVariants(q, s, PromoterVariants(0, 0))
    checkTrue(length(current) == 0L)
}

