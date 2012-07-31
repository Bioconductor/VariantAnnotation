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
    alt <- DNAStringSet(c("G", "T", "A"))
    current <- fun(query, cdsbytx, Hsapiens, alt)
    expected <- GRanges()
    checkIdentical(current, expected)
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
    current <- suppressWarnings(fun(query, cdsbytx[1:2], Hsapiens, variant))

    current_varaa <- values(current[names(current) == "B"])[["VARAA"]]
    checkTrue(as.character(current_varaa) == "")

    current_consequence <- values(current[names(current) == "B"])[["CONSEQUENCE"]]
    checkTrue(current_consequence == "not translated")

    ## TODO : add test for codon width based on 1,2,3 position
}

test_refLocsToLocalLocs <- function()
{
    ## both in 'first' cds
    query <- GRanges(seqnames="chr1",
              ranges=IRanges(rep(c(10002, 10005), 2), width=1),
              strand=c("+", "+", "-", "-"))
    current <- refLocsToLocalLocs(query, cdsbytx=cdsbytx[c(1,3)])
    expected <- IRanges(c(2, 5, 9, 6), width=1) 
    checkIdentical(values(current)[["CDSLOC"]], expected)
    expected <- IntegerList(1, 2, 3, 2) 
    checkIdentical(values(current)[["PROTEINLOC"]], expected)

    ## one in each cds
    query <- GRanges(seqnames="chr1",
                     ranges=IRanges(rep(c(10002, 10011), 2), width=1),
                     strand=c("+", "+", "-", "-"))
    current <- refLocsToLocalLocs(query, cdsbytx=cdsbytx[c(1,3)])
    expected <- IRanges(c(2, 7, 9, 4), width=1) 
    checkIdentical(values(current)[["CDSLOC"]], expected)
    expected <- IntegerList(1, 3, 3, 2) 
    checkIdentical(values(current)[["PROTEINLOC"]], expected)

    ## both in 'last' cds
    query <- GRanges(seqnames="chr1",
                     ranges=IRanges(rep(c(10010, 10013), 2), width=1),
                     strand=c("+", "+", "-", "-"))
    current <- refLocsToLocalLocs(query, cdsbytx=cdsbytx[c(1,3)])
    expected <- IRanges(c(6, 9, 5, 2), width=1) 
    checkIdentical(values(current)[["CDSLOC"]], expected)
    expected <- IntegerList(2, 3, 2, 1) 
    checkIdentical(values(current)[["PROTEINLOC"]], expected)
} 

test_predictCoding_strand <- function()
{
    variant=DNAStringSet(c("G", "G", "C", "T", "G"))
    query <- GRanges(seqnames="chr1",
              ranges=IRanges(c(rep(10003, 3), 10011, 10101), width=1),
              strand=c("+", "-", "*", "*", "*"),
              variant=variant)
    names(query) <- LETTERS[1:5]
    current <- suppressWarnings(fun(query, cdsbytx, Hsapiens, variant))

    expected <- c("G", "C", "C", "C", "G", "G", "T", "A", "C")
    checkIdentical(as.character(values(current)[["varAllele"]]), expected)

    ## ignore.strand
    strand(query) <- "+"
    p1 <- suppressWarnings(fun(query, cdsbytx, Hsapiens, variant,
        ignore.strand=TRUE))
    checkIdentical(12L, length(p1))
    p2 <- suppressWarnings(fun(query, cdsbytx, Hsapiens, variant,
        ignore.strand=FALSE))
    checkIdentical(4L, length(p2))
}



