### =========================================================================
### CollapsedVCF class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###
### See VCF(..., collapsed=TRUE) in methods-VCF-class.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters

### alt() setter
setReplaceMethod("alt", c("CollapsedVCF", "DNAStringSetList"),
    function(x, value)
{
    if (length(value) != length(rowData(x)))
        stop("length(value) must equal length(rowData(x))")
    slot(x, "fixed")$ALT <- value
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isSNV family 
###

.logicalListToVector <- function(res, alt, singleAltOnly) { 
    lst <- relist(res, alt)
    if (singleAltOnly)
        all(lst) & elementLengths(lst) == 1
    else
        any(lst)
}

.dispatchSNV_CollapsedVCF <- function(FUN, x, singleAltOnly)
{
    if (is(alt <- alt(x), "CharacterList"))
        stop("'alt' must be non-structural nucleotide values")
    res <- FUN(rep(ref(x), elementLengths(alt)),
               unlist(alt, use.names=FALSE))
    .logicalListToVector(res, alt, singleAltOnly)
}
 
setMethod("isSNV", "CollapsedVCF", 
    function(x, ..., singleAltOnly=TRUE) 
        .dispatchSNV_CollapsedVCF(.isSNV, x, singleAltOnly) 
)

setMethod("isInsertion", "CollapsedVCF", 
    function(x, ..., singleAltOnly=TRUE) 
        .dispatchSNV_CollapsedVCF(.isInsertion, x, singleAltOnly) 
)

setMethod("isDeletion", "CollapsedVCF", 
    function(x, ..., singleAltOnly=TRUE) 
        .dispatchSNV_CollapsedVCF(.isDeletion, x, singleAltOnly) 
)

setMethod("isIndel", "CollapsedVCF", 
    function(x, ..., singleAltOnly=TRUE) 
        .dispatchSNV_CollapsedVCF(.isIndel, x, singleAltOnly) 
)

setMethod("isTransition", "CollapsedVCF", 
    function(x, ..., singleAltOnly=TRUE) 
        .dispatchSNV_CollapsedVCF(.isTransition, x, singleAltOnly) 
)

setMethod("isSubstitution", "CollapsedVCF", 
    function(x, ..., singleAltOnly=TRUE) 
        .dispatchSNV_CollapsedVCF(.isSubstitution, x, singleAltOnly) 
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(show, "CollapsedVCF",
    function(object)
{
    .showVCFSubclass(object)
})
