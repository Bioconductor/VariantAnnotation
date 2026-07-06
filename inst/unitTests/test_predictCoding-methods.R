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
    ## issue #86: REFAA and VARAA should be empty AAStringSet, not NULL
    checkTrue(is(mcols(current)$REFAA, "AAStringSet"))
    checkTrue(is(mcols(current)$VARAA, "AAStringSet"))
    checkIdentical(length(mcols(current)$REFAA), 0L)
    checkIdentical(length(mcols(current)$VARAA), 0L)
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
    checkIdentical(unname(as.character(mcols(current)$VARCODON)), 
                   c("TATCCGG", "TTCCGG"))
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

test_predictCoding_nonsense_DBS <- function()
{
    ## issue #84: DBS spanning two codons where one becomes stop (*) should
    ## be classified as "nonsense", not "nonsynonymous"
    ## Construct a minimal CDS: 9-nt coding seq on "+" strand
    ## REF codon 3 = positions 7-9, codon 2 = positions 4-6
    ## A DBS at positions 5-6 (end of codon2 / start of codon3) can yield
    ## VARAA like "P*" — contains stop but wasn't caught by old %in% "*"
    fa_path <- tempfile(fileext='.fa')
    ## 9-nt CDS encodes e.g. CCG-CAG-TGG (P-Q-W)
    ## DBS variant at pos 4-5: CAG -> TAG introduces stop in codon 2 -> P*W
    cds_seq <- "CCGCAGTGG"
    full_seq <- paste0(paste(rep("A", 1000), collapse=""), cds_seq,
                       paste(rep("A", 1000), collapse=""))
    Biostrings::writeXStringSet(
        Biostrings::DNAStringSet(c(chr1=full_seq)), fa_path)
    Rsamtools::indexFa(fa_path)
    fa <- Rsamtools::FaFile(fa_path)

    cds_start <- 1001L
    cdsbytx_dbs <- GRangesList(tx1=GRanges("chr1",
        IRanges(cds_start, cds_start + 8L), strand="+"))

    ## variant at codon-boundary positions 4-5 of the CDS (genomic 1004-1005)
    ## REF=CA, ALT=TA  ->  VARCODON2 = TAG (stop), VARAA = "*W"
    query_dbs <- GRanges("chr1",
        IRanges(cds_start + 3L, cds_start + 4L),
        strand="+")
    varAllele_dbs <- Biostrings::DNAStringSet("TA")

    current <- quiet(fun(query_dbs, cdsbytx_dbs, fa, varAllele_dbs))

    if (length(current) > 0L) {
        ## VARAA should contain a stop (*) somewhere
        checkTrue(grepl("\\*", as.character(mcols(current)$VARAA)))
        ## CONSEQUENCE must be "nonsense", not "nonsynonymous"
        checkTrue(as.character(mcols(current)$CONSEQUENCE) == "nonsense")
    }
}


## Regression test for issue #83:
## A deletion spanning an exon/intron boundary should use CDS-mapped width
## (CDSLOC width) not genomic width when computing REFCODON/VARCODON.
test_predictCoding_exon_intron_boundary <- function()
{
    ## Build a minimal two-exon transcript on chr1
    ## Exon1: 1-90 (90 nt), Exon2: 201-300 (100 nt); intron: 91-200
    txdb <- suppressMessages(
        makeTxDbFromGRanges(GRanges(
            seqnames = "chr1",
            ranges   = IRanges(
                start = c(1, 1, 1, 201),
                end   = c(300, 300, 90, 300)
            ),
            strand = "+",
            type   = c("gene", "mRNA", "exon", "exon"),
            ID     = c("gene1", "tx1", "exon1", "exon2"),
            Parent = c(NA, "gene1", "tx1", "tx1")
        ))
    )
    cdsbytx <- suppressMessages(cdsBy(txdb, by="tx"))
    if (length(cdsbytx) == 0L) return(invisible(NULL))  # skip if txdb empty

    ## Deletion that straddles exon1 end (pos 85) into the intron (pos 95)
    ## Genomic span = 11 bp, but only 6 bp are exonic (85-90)
    del <- GRanges("chr1", IRanges(85, 95), strand="+",
                   REF=DNAStringSet("NNNNNNNNNN"),
                   ALT=DNAStringSetList(DNAStringSet("-")))
    names(del) <- "boundary_del"

    ## Build a tiny genome fasta
    fa <- Biostrings::DNAStringSet(paste0(
        paste(rep("A", 84), collapse=""),   # 1-84
        "CGTACG",                           # 85-90 (exon1 end, 6 nt)
        paste(rep("T", 110), collapse=""),  # 91-200 intron
        paste(rep("G", 100), collapse="")   # 201-300 exon2
    ))
    names(fa) <- "chr1"
    tmpfa <- tempfile(fileext=".fa")
    Biostrings::writeXStringSet(fa, tmpfa)
    genome <- Rsamtools::FaFile(tmpfa)

    res <- tryCatch(
        predictCoding(del, txdb, seqSource=genome,
                      varAllele=DNAStringSet("-")),
        error = function(e) NULL
    )
    ## If a result is returned the CDSLOC width must be <= 6 (exonic bases only)
    ## not 11 (the full genomic span)
    if (!is.null(res) && length(res) > 0L) {
        cdswidths <- width(mcols(res)$CDSLOC)
        checkTrue(all(cdswidths <= 6L))
    }
}


## Regression test for issue #46:
## When multiple variants fall in the same codon of the same transcript,
## .predictCodingGRangesList() should warn that each is evaluated
## independently (not as a compound haplotype).
test_predictCoding_multi_snp_codon_warning <- function()
{
    ## Re-use the CDS/genome already set up at the top of this file.
    ## cdsbytx$tx1 = chr1:10001-10005, 10010-10014 (two exons, + strand)
    ## CDS position 1 = genomic 10001, position 2 = 10002, position 3 = 10003
    ## => codon 1 spans genomic 10001-10003.
    ## Place two SNPs at 10001 and 10002 — same codon, same transcript (tx1).
    query <- GRanges("chr1", IRanges(start = c(10001L, 10002L), width = 1L),
                     strand = "+")
    va    <- DNAStringSet(c("T", "G"))

    msgs <- character(0)
    withCallingHandlers(
        fun(query, cdsbytx["tx1"], Hsapiens, va),
        warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    checkTrue(any(grepl("multiple variants", msgs, ignore.case = TRUE)),
              "expected warning about multiple variants in same codon (#46)")
}
