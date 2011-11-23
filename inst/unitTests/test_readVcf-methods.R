f1 <- system.file("extdata", "ex1.vcf", package="VariantAnnotation")
f2 <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")

test_ranges <- function()
{
    vcf <- readVcf(f1, "hg18")
    checkEquals(width(rowData(vcf)), width(values(rowData(vcf))[["REF"]]))
    checkEquals(scanVcf(f1)[[1]]$POS, start(rowData(vcf))) 
}

test_metadata <- function()
{
    vcf <- readVcf(f1, "hg18")
    checkTrue(all(rownames(exptData(vcf)[["HEADER"]][["FORMAT"]]) %in% 
              names(assays(vcf)))) 
    checkTrue(all(rownames(exptData(vcf)[["HEADER"]][["INFO"]]) %in% 
              names(values(rowData(vcf)))))
    checkIdentical("hg18", unique(genome(rowData(vcf)))) 
}

test_accessors <- function()
{
    vcf <- readVcf(f1, "hg18")
    checkTrue(class(vcf) == "SummarizedExperiment")
    checkTrue(class(rowData(vcf)) == "GRanges")
    checkTrue(class(assays(vcf)) == "SimpleList")
    checkTrue(class(exptData(vcf)) == "SimpleList")
    checkTrue(class(colData(vcf)) == "DataFrame")
}

