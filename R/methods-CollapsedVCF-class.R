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

.logicalListToVector <- function(res, alt, singleAltOnly, ignore) { 
    lst <- relist(res, alt)
    if (singleAltOnly && !any(ignore))
        all(lst) & elementLengths(lst) == 1
    else
        any(lst)
}

.dispatchSNV_CollapsedVCF <- function(FUN, x, singleAltOnly)
{
    flat <- unlist(alt(x), use.names=FALSE)
    ignore <- logical(length(flat))
    if (is(alt(x), "CharacterList")) {
        ## handle gvcf, not structural
        ignore <- grepl("NON_REF", flat, fixed=TRUE)
        if (!any(ignore) && any(.isStructural(flat)))
            stop("'alt' values must be nucleotides")
    } 
    res <- FUN(rep(ref(x), elementLengths(alt(x))), flat, ignore)
    .logicalListToVector(res, alt(x), singleAltOnly, ignore)
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
