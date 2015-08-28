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
    if (length(value) != length(rowRanges(x)))
        stop("length(value) must equal length(rowRanges(x))")
    slot(x, "fixed")$ALT <- value
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isSNV family 
###

.dispatchSNV_CollapsedVCF <- function(FUN, x, singleAltOnly)
{
    alt <- alt(x)
    flat <- unlist(alt, use.names=FALSE)
    gvcf <- structural <- NULL
    if (is(flat, "character")) {
        gvcf <- grepl("NON_REF", flat, fixed=TRUE)
        structural <- .isStructural(flat)
    }
    res <- FUN(rep(ref(x), elementLengths(alt(x))), flat)
    if (!is.null(structural))
        res[structural | is.na(structural)] <- FALSE
    if (any(gvcf))
        res[gvcf] <- TRUE
    ## relist
    lst <- relist(res, alt)
    if (singleAltOnly)
        all(lst) & elementLengths(lst) == 1
    else
        any(lst)
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

setMethod("isDelins", "CollapsedVCF", 
    function(x, ..., singleAltOnly=TRUE) 
        .dispatchSNV_CollapsedVCF(.isDelins, x, singleAltOnly) 
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
