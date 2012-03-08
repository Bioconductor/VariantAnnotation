setMethod(ScanVcfParam, "ANY",
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, ...)
{
    ScanBcfParam(fixed, info, geno, trimEmpty, which=which,
                 class="ScanVcfParam")
})

setMethod(ScanVcfParam, "missing",
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, ...)
{
    ScanBcfParam(fixed, info, geno, trimEmpty, class="ScanVcfParam")
})

## accessors

vcfFixed <- function(object) slot(object, "fixed")

vcfInfo <- function(object) slot(object, "info")

vcfGeno <- function(object) slot(object, "geno")

vcfTrimEmpty <- function(object) slot(object, "trimEmpty")

vcfWhich <- function(object) 
{
    rl <- slot(object, "which")
    gr <- as(rl, "GRanges")
    names(gr) <- names(unlist(rl, use.names=FALSE))
    gr
}
