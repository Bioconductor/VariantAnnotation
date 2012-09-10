grl <- GRangesList(
    A1=GRanges("20", IRanges(14370, width=1)),
    B=GRanges("20", IRanges(17330, width=1)),
    A2=GRanges("20", IRanges(14370, width=2)),
    C1=GRanges("20", IRanges(1110696, width=1)), 
    C2=GRanges("20", IRanges(1110696, width=2)), 
    C3=GRanges("20", IRanges(1110696, width=3))) 

test_summarizeOverlaps_findOverlaps <- function()
{
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
    vcf <- readVcf(fl, "hg19")
    sv <- summarizeVariants(grl, vcf, findOverlaps)
    current <- unname(colSums(assays(sv)$counts))
    target <- c(3, 6, 5)
    checkIdentical(current, target) 

    genes <- factor(c("A", "B", "A", "C", "C", "C"))
    sv <- summarizeVariants(grl, vcf, findOverlaps, subjectFactor=genes)
    current <- unname(colSums(assays(sv)$counts))
    target <- c(1, 3, 2)
    checkIdentical(current, target) 
    checkTrue(nrow(assays(sv)$counts) == length(levels(genes)))
}

