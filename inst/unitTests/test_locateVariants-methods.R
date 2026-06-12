library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
cdsbytx <- cdsBy(txdb, use.names=TRUE)
intbytx <- intronsByTranscript(txdb)
txbygene <- transcriptsBy(txdb, "gene")

#
# As of Oct 16 2025, 3.22 preparation, we have the following
# but the test positions in gr below for locateVariants are
# inconsistent with these addresses
#
#> cdsbytx[1:3]
#GRangesList object of length 3:
#$ENST00000641515.2_6
#GRanges object with 2 ranges and 3 metadata columns:
#      seqnames      ranges strand |    cds_id    cds_name exon_rank
#         <Rle>   <IRanges>  <Rle> | <integer> <character> <integer>
#  [1]     chr1 65565-65573      + |         1        <NA>         2
#  [2]     chr1 69037-70008      + |         2        <NA>         3
#  -------
#  seqinfo: 298 sequences (2 circular) from hg19 genome
#
#$ENST00000426406.1
#GRanges object with 1 range and 3 metadata columns:
#      seqnames        ranges strand |    cds_id    cds_name exon_rank
#         <Rle>     <IRanges>  <Rle> | <integer> <character> <integer>
#  [1]     chr1 367659-368597      + |         3        <NA>         1
#  -------
#  seqinfo: 298 sequences (2 circular) from hg19 genome
#
#$ENST00000616016.5_7
#GRanges object with 14 ranges and 3 metadata columns:
#       seqnames        ranges strand |    cds_id    cds_name exon_rank
#          <Rle>     <IRanges>  <Rle> | <integer> <character> <integer>
#   [1]     chr1 859812-860328      + |         4        <NA>         1
#   [2]     chr1 861302-861393      + |         5        <NA>         2
#   [3]     chr1 865535-865716      + |         7        <NA>         3
#   [4]     chr1 866419-866469      + |         9        <NA>         4
#   [5]     chr1 871152-871276      + |        11        <NA>         5
#

gr <- GRanges("chr22", 
    IRanges(c(16268137, 16287254, 16190792, 16164570,
              18209442, 18121652, 24314750, 25508661), 
            width=c(1,1,1,1,3,3,2,2)),
    strand=c("-", "-", "-", "+", "+", "+", "+", "+"))

test_locateVariants_upstream_downstream <- function()
{
#    loc <- locateVariants(gr, txdb, IntergenicVariants(1, 1))  # VJC 10/16/2025 not true
#    target <- CharacterList(character(), character())
#    checkIdentical(loc$FOLLOWID, target)

    loc <- locateVariants(gr, txbygene, IntergenicVariants(100, 100))
#    target <- CharacterList(character(), "100037417")
    target <- CharacterList("100037417")
    checkIdentical(loc$FOLLOWID, target) 

    loc <- locateVariants(gr, txbygene, IntergenicVariants(100000, 100000))
# 100037417 100652871 4282 66035
    target <- CharacterList(c("100037417", "100652871", "4282", "66035"))
    checkIdentical(loc$FOLLOWID, target)
#  23523,2953
    target <- CharacterList(c("23523", "2953")) # , "391322"))
    checkIdentical(loc$PRECEDEID, target)
}

test_locateVariants_queryAsVCF <- function()
{
    fl <- system.file("extdata", "gl_chr1.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    seqlevels(vcf) <- c("1" = "chr1")
    loc1 <- locateVariants(vcf, txdb, IntergenicVariants())
    loc2 <- locateVariants(rowRanges(vcf), txdb, IntergenicVariants())
    checkIdentical(loc1, loc2) 
}

test_locateVariants_ignore.strand <- function()
{
    cdsbytx <- cdsbytx[1:5]
#    gr <- GRanges("chr1", IRanges(c(12190, 12595, 13403), width=1), "-")
    gr = GRanges("chr1", IRanges(c(65565, 367659, 859812), width=1), "-")
    loc1 <- locateVariants(gr, cdsbytx, CodingVariants(), 
                           ignore.strand=TRUE)
    checkIdentical(c(1L, 2L, 3L, 3L), mcols(loc1)$QUERYID)   # changed VJC 16 oct 2025
    loc2 <- locateVariants(gr, cdsbytx, CodingVariants(), 
                           ignore.strand=FALSE)
    checkIdentical(integer(), mcols(loc2)$QUERYID) 
    loc1 <- locateVariants(gr, cdsbytx, SpliceSiteVariants(), 
                           ignore.strand=TRUE)
print(loc1)
    checkIdentical(c(1L, 2L, 3L, 3L), mcols(loc1)$QUERYID) 
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

## ----------------------------------------------------------------------------
## Issue #81: large deletions extending beyond transcript bounds were silently
## dropped by locateVariants() / predictCoding() because mapToTranscripts()
## uses findOverlaps(type="within") internally.
## ----------------------------------------------------------------------------
test_locateVariants_large_deletion_not_dropped <- function()
{
    ## Build a synthetic 3-exon transcript on a 2000 bp chromosome
    ## Exon1: 100-299, Exon2: 400-599, Exon3: 700-899  (all "+")
    cds_ranges <- GRanges("chr1",
        IRanges(start=c(100L, 400L, 700L),
                end  =c(299L, 599L, 899L)),
        strand="+")
    subject <- GRangesList(tx1=cds_ranges)

    ## Small variant fully within exon2 -> should be found
    q_small <- GRanges("chr1", IRanges(450L, 460L), strand="+")
    loc_small <- locateVariants(q_small, subject, CodingVariants())
    checkTrue(length(loc_small) > 0L,
        "small intra-exon variant must be located (sanity check)")

    ## Large deletion spanning ALL of exon2 and into the flanking introns:
    ## genomic 350-650, which extends beyond the exon2 bounds (400-599)
    ## Before the fix this was silently dropped.
    q_large <- GRanges("chr1", IRanges(350L, 650L), strand="+")
    loc_large <- locateVariants(q_large, subject, CodingVariants())
    checkTrue(length(loc_large) > 0L,
        "large deletion spanning exon bounds must not be silently dropped (#81)")

    ## QUERYID should refer back to query index 1
    checkIdentical(unique(loc_large$QUERYID), 1L)
}

test_localCoordinates_large_deletion_not_dropped <- function()
{
    ## Issue #81: test via .localCoordinates() (the internal function called by
    ## predictCoding) that large deletions spanning CDS bounds are not dropped.
    cds_start <- 501L
    cds_ranges <- GRanges("chr1",
        IRanges(cds_start, cds_start + 8L), strand="+",
        cds_id=1L, cds_name=NA_character_, exon_rank=1L)
    cdsbytx <- GRangesList(tx1=cds_ranges)

    q_large <- GRanges("chr1",
        IRanges(cds_start - 10L, cds_start + 18L), strand="+")
    mcols(q_large)$varAllele <- Biostrings::DNAStringSet("A")

    map <- VariantAnnotation:::.localCoordinates(
        q_large, cdsbytx, ignore.strand=FALSE)
    checkTrue(length(map) > 0L,
        ".localCoordinates must not drop large deletion spanning CDS bounds (#81)")
}
