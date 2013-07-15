### =========================================================================
### VRangesList: Where there is a Ranges, there must be a List
### -------------------------------------------------------------------------
###

VRangesList <- function(x) {
  new("SimpleVRangesList", GenomicRangesList(x))
}

setMethod("stackSamples", "SimpleVRangesList", function(x) {
  stacked <- unlist(x, use.names=FALSE)
  if (!is.null(names(x)))
    sampleNames(stacked) <- Rle(names(x), elementLengths(x))
  stacked
})

setMethod("merge", c("SimpleVRangesList", "missing"), function(x, y) {
  ans <- unique(unlist(x, use.names=FALSE))
  ad.sum <- 0L
  td.sum <- 0L
  rd.sum <- 0L
  for (xi in as.list(x, use.names=FALSE)) {
    m <- match(ans, xi)
    td.sum <- td.sum + ifelse(is.na(m), 0L, totalDepth(xi)[m])
    ad.sum <- ad.sum + ifelse(is.na(m), 0L, altDepth(xi)[m])
    rd.sum <- rd.sum + ifelse(is.na(m), 0L, refDepth(xi)[m])
  }
  altDepth(ans) <- ad.sum
  totalDepth(ans) <- td.sum
  refDepth(ans) <- rd.sum
  ans
})
