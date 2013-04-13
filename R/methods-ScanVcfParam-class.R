### =========================================================================
### ScanVcfParam class methods 
### =========================================================================


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Constructor 
##

setMethod(ScanVcfParam, "ANY",
    function(fixed=character(), info=character(), geno=character(), 
             samples=character(), trimEmpty=TRUE, which, ...)
{
    ScanBcfParam(fixed, info, geno, samples, trimEmpty=trimEmpty, 
                 which=which, class="ScanVcfParam")
})

setMethod(ScanVcfParam, "missing",
    function(fixed=character(), info=character(), geno=character(), 
             samples=character(), trimEmpty=TRUE, which, ...)
{
    ScanBcfParam(fixed, info, geno, samples, trimEmpty=trimEmpty,
                 class="ScanVcfParam")
})

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Validity 
##

.valid.ScanVcfParam <- function(x)
{
    samples <- vcfSamples(x)
    geno <- vcfGeno(x)
    if (any(is.na(samples)) && length(geno) > 0L)
        return("ScanVcfParam: 'geno' cannot be specified if 'samples' is 'NA'")
    if (any(is.na(geno)) && length(samples) > 0L)
        return("ScanVcfParam: 'samples' cannot be specified if 'geno' is 'NA'")

    NULL 
}

setValidity("ScanVcfParam", .valid.ScanVcfParam,
    where=asNamespace("VariantAnnotation"))

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

vcfSamples <- function(object) slot(object, "samples")
"vcfSamples<-" <- function(object, value) 
{
    slot(object, "samples") <- value
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
