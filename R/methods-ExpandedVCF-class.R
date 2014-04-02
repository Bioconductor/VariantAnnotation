### =========================================================================
### ExpandedVCF class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###
### See VCF(..., collapsed=FALSE) in methods-VCF-class.R
 
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

### alt() setter
setReplaceMethod("alt", c("ExpandedVCF", "DNAStringSet"),
    function(x, value)
{
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

VRangesForMatching <- function(x) {
    with(rowData(x), VRanges(seqnames, IRanges(start, end), ref=REF, alt=ALT))
}

setMethod("match", c("ExpandedVCF", "ExpandedVCF"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL,
             method = c("auto", "quick", "hash"))
{
    x <- VRangesForMatching(x)
    table <- VRangesForMatching(table)
    callGeneric()
})

setMethod("match", c("ExpandedVCF", "VRanges"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL,
             method = c("auto", "quick", "hash"))
{
    x <- VRangesForMatching(x)
    callGeneric()
})

setMethod("match", c("VRanges", "ExpandedVCF"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL,
             method = c("auto", "quick", "hash"))
{
    table <- VRangesForMatching(table)
    callGeneric()
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(show, "ExpandedVCF",
    function(object)
{
    .showVCFSubclass(object)
})
