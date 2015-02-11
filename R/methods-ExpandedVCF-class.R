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
### isSNV family 
###

.dispatchSNV_ExpandedVCF <- function(FUN, x)
{
    alt <- alt(x)
    ignore <- logical(length(alt))
    if (is(alt, "character")) {
        ## handle gvcf, not structural
        ignore <- grepl("NON_REF", unlist(alt, use.names=FALSE), fixed=TRUE)
        if (!any(ignore))
            if (any(.isStructural(alt)))
                stop("'alt' values must be nucleotides")
        if (any(elementLengths(alt) > 1L))
            stop("all elementLengths(alt) must be 1")
    }
    FUN(ref(x), alt, ignore)
}

setMethod("isSNV", "ExpandedVCF", 
    function(x, ...) 
        .dispatchSNV_ExpandedVCF(.isSNV, x) 
)

setMethod("isInsertion", "ExpandedVCF", 
    function(x, ...) 
        .dispatchSNV_ExpandedVCF(.isInsertion, x) 
)

setMethod("isDeletion", "ExpandedVCF", 
    function(x, ...) 
        .dispatchSNV_ExpandedVCF(.isDeletion, x) 
)

setMethod("isIndel", "ExpandedVCF", 
    function(x, ...) 
        .dispatchSNV_ExpandedVCF(.isIndel, x) 
)

setMethod("isDelins", "ExpandedVCF", 
    function(x, ...) 
        .dispatchSNV_ExpandedVCF(.isDelins, x) 
)

setMethod("isTransition", "ExpandedVCF",
    function(x, ...) 
        .dispatchSNV_ExpandedVCF(.isTransition, x) 
)

setMethod("isSubstitution", "ExpandedVCF",
    function(x, ...) 
        .dispatchSNV_ExpandedVCF(.isSubstitution, x) 
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(show, "ExpandedVCF",
    function(object)
{
    .showVCFSubclass(object)
})
