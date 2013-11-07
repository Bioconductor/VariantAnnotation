### =========================================================================
### CollapsedVCF class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###

CollapsedVCF <- 
    function(rowData=GRanges(), colData=DataFrame(), 
             exptData=SimpleList(header=VCFHeader()), 
             fixed=DataFrame(),
             info=DataFrame(row.names=seq_along(rowData)), 
             geno=SimpleList(),
             ..., verbose=FALSE)
{
    VCF(rowData, colData, exptData, fixed, info, geno,
        collapsed=TRUE, verbose=verbose)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

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
### show
###

setMethod(show, "CollapsedVCF",
    function(object)
{
    .showVCFSubclass(object)
})
