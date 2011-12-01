f1 <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
f2 <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
f3 <- system.file("extdata", "structural.vcf", package="VariantAnnotation")

test_readVcf_ranges <- function()
{
    vcf <- readVcf(f1, "hg19")
    checkEquals(width(rowData(vcf)), width(values(rowData(vcf))[["REF"]]))
    checkEquals(scanVcf(f1)[[1]]$POS, start(rowData(vcf))) 
}

test_readVcf_metadata <- function()
{
    vcf <- readVcf(f1, "hg19")
    checkTrue(all(rownames(exptData(vcf)[["HEADER"]][["FORMAT"]]) %in% 
              names(assays(vcf)))) 
    checkIdentical("hg19", unique(genome(rowData(vcf)))) 
}

test_readVcf_accessors <- function()
{
    vcf <- readVcf(f1, "hg19")
    checkTrue(class(vcf) == "SummarizedExperiment")
    checkTrue(class(rowData(vcf)) == "GRanges")
    checkTrue(class(assays(vcf)) == "SimpleList")
    checkTrue(class(exptData(vcf)) == "SimpleList")
    checkTrue(class(colData(vcf)) == "DataFrame")

    vcf_st <- readVcf(f3, "hg19")
    checkTrue(class(rowData(vcf_st)) == "GRanges")
}

