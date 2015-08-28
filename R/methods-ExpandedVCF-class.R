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
    if (length(value) != length(rowRanges(x)))
        stop("length(value) must equal length(rowRanges(x))")
    slot(x, "fixed")$ALT <- value
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

VRangesForMatching <- function(x) {
    with(rowRanges(x), VRanges(seqnames, IRanges(start, end), ref=REF, alt=ALT))
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
    flat <- unlist(alt, use.names=FALSE)
    gvcf <- structural <- NULL
    if (is(flat, "character")) {
        gvcf <- grepl("NON_REF", flat, fixed=TRUE)
        structural <- .isStructural(flat)
    }
    res <- FUN(ref(x), alt)
    if (!is.null(structural))
        res[structural | is.na(structural)] <- FALSE
    if (any(gvcf))
        res[gvcf] <- TRUE

    res
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
