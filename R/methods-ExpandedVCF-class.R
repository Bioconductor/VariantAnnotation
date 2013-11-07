### =========================================================================
### ExpandedVCF class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###

ExpandedVCF <- 
    function(rowData=GRanges(), colData=DataFrame(), 
             exptData=SimpleList(header=VCFHeader()), 
             fixed=DataFrame(),
             info=DataFrame(row.names=seq_along(rowData)), 
             geno=SimpleList(),
             ..., verbose=FALSE)
{
    VCF(rowData, colData, exptData, fixed, info, geno,
        collapsed=FALSE, verbose=verbose)
}
 
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
### show
###

setMethod(show, "ExpandedVCF",
    function(object)
{
    .showVCFSubclass(object)
})
