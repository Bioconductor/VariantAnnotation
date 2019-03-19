library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
cdsbytx <- cdsBy(txdb, use.names=TRUE)
intbytx <- intronsByTranscript(txdb)
txbygene <- transcriptsBy(txdb, "gene")

gr <- GRanges("chr22", 
    IRanges(c(16268137, 16287254, 16190792, 16164570,
              18209442, 18121652, 24314750, 25508661), 
            width=c(1,1,1,1,3,3,2,2)),
    strand=c("-", "-", "-", "+", "+", "+", "+", "+"))

test_locateVariants_upstream_downstream <- function()
{
    loc <- locateVariants(gr, txdb, IntergenicVariants(1, 1))
    target <- CharacterList(character(), character())
    checkIdentical(loc$FOLLOWID, target)

    loc <- locateVariants(gr, txbygene, IntergenicVariants(2, 2))
    target <- CharacterList(character(), "100037417")
    checkIdentical(loc$FOLLOWID, target) 

    loc <- locateVariants(gr, txbygene, IntergenicVariants(100000, 100000))
    target <- CharacterList("23784", c("100037417","4282", "66035"))
    checkIdentical(loc$FOLLOWID, target)
    target <- CharacterList(character(), c("23523", "2953", "391322"))
    checkIdentical(loc$PRECEDEID, target)
}

test_locateVariants_queryAsVCF <- function()
{
    fl <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    vcf <- renameSeqlevels(vcf, c("1" = "chr1"))
    loc1 <- locateVariants(vcf, txdb, IntergenicVariants())
    loc2 <- locateVariants(rowRanges(vcf), txdb, IntergenicVariants())
    checkIdentical(loc1, loc2) 
}

test_locateVariants_ignore.strand <- function()
{
    cdsbytx <- cdsbytx[1:5]
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

    q <- GRanges("chr22", IRanges(50310410, 50310420))
    current <- locateVariants(q, txdb, PromoterVariants())
    checkIdentical(unique(current$GENEID), "79174")
}

test_locateVariants_match_predictCoding <- function()
{
    library(BSgenome.Hsapiens.UCSC.hg19)
    gr <- GRanges("chr20", IRanges(
        start=c(77055, 77054, 77054, 77058, 77057, 77057, 77055), 
        end=c(77055, 77055, 77055, 77058, 77058, 77058, 77054)),
        paramRangeID=rep(NA, 7))
    fixed <- DataFrame(
        REF=DNAStringSet(c('T', 'AT', 'AT', 'A', 'AA', 'AA', 'T')), 
        ALT=DNAStringSetList('G', 'A', 'ATT', 'G', 'A', 'AAT', 'G'), 
        QUAL=70, FILTER="PASS")
    vcf <- VCF(rowRanges=gr, fixed=fixed)

    ## coding regions match, zero-width
    loc <- locateVariants(vcf, txdb, CodingVariants())
    coding <- predictCoding(vcf, txdb, Hsapiens)
    checkIdentical(loc$QUERYID, as.integer(1:7))
    checkIdentical(length(coding) , 6L)
    checkIdentical(unname(loc$CDSID[1:6]), unname(coding$CDSID))
    checkIdentical(unname(as.character(coding$VARCODON[c(1,4)])), 
                   as.character(DNAStringSet(c("AAG", "TAG"))))
}
    gr2 <- GRanges("chr20", IRanges(
        start=c(5, 77054, 77054, 77058, 77057, 77057, 77055), 
        end=c(55, 77055, 77055, 77058, 77058, 77058, 77054)),
        paramRangeID=rep(NA, 7))
