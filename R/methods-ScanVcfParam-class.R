setMethod(ScanVcfParam, "ANY",
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, asGRanges=FALSE, ...)
{
    ScanBcfParam(fixed, info, geno, trimEmpty, which, asGRanges, 
                 class="ScanVcfParam")
})

setMethod(ScanVcfParam, "missing",
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, asGRanges=FALSE, ...)
{
    ScanBcfParam(fixed, info, geno, trimEmpty, asGRanges=asGRanges, 
                 class="ScanVcfParam")
})

## accessors

vcfFixed <- function(object) slot(object, "fixed")

vcfInfo <- function(object) slot(object, "info")

vcfGeno <- function(object) slot(object, "geno")

vcfTrimEmpty <- function(object) slot(object, "trimEmpty")

vcfWhich <- function(object) as(slot(object, "which"), "GRanges")

vcfAsGRanges <- function(object) slot(object, "asGRanges")

