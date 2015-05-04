quiet <- suppressWarnings
library(BSgenome.Hsapiens.UCSC.hg19)
fun <- VariantAnnotation:::.predictCodingGRangesList
cdsbytx <- GRangesList(tx1=GRanges(seqnames="chr1", 
                                   IRanges(c(10001, 10010), width=5), 
                                   strand="+"),
                       tx2=GRanges(seqnames="chr1", 
                                   IRanges(c(10100, 10001), width=5), 
                                   strand="-"),
                       tx3=GRanges(seqnames="chr1", 
                                   IRanges(c(10010, 10001), width=5), 
                                   strand="-"))

test_predictCoding_empty <- function()
{
    query <- GRanges("chr1", IRanges(start=c(1, 10, 20), width=1))
    current <- fun(query, cdsbytx, Hsapiens, DNAStringSet(c("G", "T", "A")))
    checkIdentical(dim(mcols(current)), c(0L, 8L))
}

test_predictCoding_varAllele <- function()
{
    variant=DNAStringSet(c("G", "", "C", "AA", "GGA"))
    query <- GRanges(seqnames="chr1",
              ranges=IRanges(c(rep(10003, 3), 10011, 10101), 
                             width=c(1, 1, 1, 2, 3)),
              strand=c("+", "-", "*", "*", "*"),
              variant=variant)
    names(query) <- LETTERS[1:5]
    current <- quiet(fun(query, cdsbytx[1:2], Hsapiens, variant))

    current_varaa <- values(current[names(current) == "B"])[["VARAA"]]
    checkTrue(as.character(current_varaa) == "")

    current_consequence <- 
      values(current[names(current) == "B"])[["CONSEQUENCE"]]
    checkTrue(current_consequence == "not translated")

    variant=DNAStringSet(c("GGA", "GGA"))
    query <- GRanges("chr1", IRanges(rep(10101, 2), width=c(2,3)), 
                     variant=variant)
    current <- quiet(fun(query, cdsbytx[1:2], Hsapiens, variant))
    checkIdentical(as.character(mcols(current)$VARCODON), c("", "TTCCGG"))
}

test_mapToTranscripts <- function()
{
    ## both in 'first' cds
    query <- GRanges(seqnames="chr1",
              ranges=IRanges(rep(c(10002, 10005), 2), width=1),
              strand=c("+", "+", "-", "-"))
    current <- mapToTranscripts(query, cdsbytx[c(1,3)], ignore.strand=FALSE)
    expected <- IRanges(c(2, 5, 9, 6), width=1) 
    checkIdentical(ranges(current), expected)

    ## one in each cds
    query <- GRanges(seqnames="chr1",
                     ranges=IRanges(rep(c(10002, 10011), 2), width=1),
                     strand=c("+", "+", "-", "-"))
    current <- mapToTranscripts(query, cdsbytx[c(1,3)], ignore.strand=FALSE)
    expected <- IRanges(c(2, 7, 9, 4), width=1) 
    checkIdentical(ranges(current), expected)

    ## both in 'last' cds
    query <- GRanges(seqnames="chr1",
                     ranges=IRanges(rep(c(10010, 10013), 2), width=1),
                     strand=c("+", "+", "-", "-"))
    current <- mapToTranscripts(query, cdsbytx[c(1,3)], ignore.strand=FALSE)
    expected <- IRanges(c(6, 9, 5, 2), width=1) 
    checkIdentical(ranges(current), expected)
} 

test_predictCoding_strand <- function()
{
    variant=DNAStringSet(c("G", "G", "C", "T", "G"))
    query <- GRanges(seqnames="chr1",
              ranges=IRanges(c(rep(10003, 3), 10011, 10101), width=1),
              strand=c("+", "-", "*", "*", "*"),
              variant=variant)
    names(query) <- LETTERS[1:5]

    current <- quiet(fun(query, cdsbytx, Hsapiens, variant))
    expected <- c("G", "C", "C", "C", "G", "G", "T", "A", "C")
    checkIdentical(as.character(mcols(current)$varAllele), expected)

    ## query "+", subject "-"
    v <- variant[2]
    q <- query[2]
    strand(q) <- "+"
    s <- cdsbytx[3]
    current <- quiet(fun(q, s,  Hsapiens, v, ignore.strand=FALSE))
    checkIdentical(length(current), 0L)

    current <- quiet(fun(q, s, Hsapiens, v, ignore.strand=TRUE))
    checkIdentical(as.character(mcols(current)$REFAA), "V")
    checkIdentical(as.character(mcols(current)$VARAA), "A")
    checkIdentical(mcols(current)$CDSLOC, IRanges(8, 8))

    ## query "-", subject "+"
    strand(q) <- "-"
    s <- cdsbytx[1]
    current <- quiet(fun(q, s, Hsapiens, v, ignore.strand=FALSE))
    checkIdentical(length(current), 0L)

    current <- quiet(fun(q, s, Hsapiens, v, ignore.strand=TRUE))
    checkIdentical(as.character(mcols(current)$REFAA), "*")
    checkIdentical(as.character(mcols(current)$VARAA), "*")
    checkIdentical(mcols(current)$CDSLOC, IRanges(3, 3))
}

