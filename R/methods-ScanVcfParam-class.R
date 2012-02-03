setMethod(ScanVcfParam, "ANY",
    function(info=character(), geno=character(), trimEmpty=TRUE, which, 
             asGRanges=FALSE, ...)
{
    ScanBcfParam(info, geno, trimEmpty, which, asGRanges, class="ScanVcfParam")
})

setMethod(ScanVcfParam, "missing",
    function(info=character(), geno=character(), trimEmpty=TRUE, which, 
             asGRanges=FALSE, ...)
{
    ScanBcfParam(info, geno, trimEmpty, asGRanges=asGRanges, class="ScanVcfParam")
})

## accessors

vcfInfo <- function(object) slot(object, "info")

vcfGeno <- function(object) slot(object, "geno")

vcfTrimEmpty <- function(object) slot(object, "trimEmpty")

vcfWhich <- function(object) as(slot(object, "which"), "GRanges")

vcfAsGRanges <- function(object) slot(object, "asGRanges")

