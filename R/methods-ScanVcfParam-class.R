### =========================================================================
### ScanVcfParam class methods 
### =========================================================================


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Constructor 
##

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

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Getters and Setters
##

vcfFixed <- function(object) slot(object, "fixed")
"vcfFixed<-" <- function(object, value) 
{
    slot(object, "fixed") <- value
    object
}

vcfInfo <- function(object) slot(object, "info")
"vcfInfo<-" <- function(object, value) 
{
    slot(object, "info") <- value
    object
}

vcfGeno <- function(object) slot(object, "geno")
"vcfGeno<-" <- function(object, value) 
{
    slot(object, "geno") <- value
    object
}

vcfTrimEmpty <- function(object) slot(object, "trimEmpty")
"vcfTrimEmpty<-" <- function(object, value) 
{
    slot(object, "trimEmpty") <- value
    object
}

vcfWhich <- function(object) 
{
    slot(object, "which")
}
"vcfWhich<-" <- function(object, value)
{
    slot(object, "which") <- as(value, "RangesList")
    object
}
