### =========================================================================
### VRangesList: Where there is a Ranges, there must be a List
### -------------------------------------------------------------------------
###

VRangesList <- function(...) {
  new("SimpleVRangesList", GenomicRangesList(...))
}

setMethod("stackSamples", "VRangesList", function(x) {
  stacked <- unlist(x, use.names=FALSE)
  if (!is.null(names(x)))
    sampleNames(stacked) <- Rle(names(x), elementNROWS(x))
  stacked
})
