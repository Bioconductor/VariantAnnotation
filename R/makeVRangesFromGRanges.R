### =========================================================================
### makeVRangesFromGRanges()
### -------------------------------------------------------------------------

makeVRangesFromGRanges <- function(gr, 
                                   ref.field="ref",
                                   alt.field="alt",
                                   totalDepth.field="totalDepth", 
                                   altDepth.field="altDepth", 
                                   refDepth.field="refDepth", 
                                   sampleNames.field="sampleNames",
                                   keep.extra.columns=TRUE)
{
    ## check args
    if (!is(gr, "GenomicRanges"))
        stop("'gr' must be a GenomicRanges object")
    args <- list(ref.field=ref.field, 
                 alt.field=alt.field, 
                 totalDepth.field=totalDepth.field, 
                 altDepth.field=altDepth.field,
                 refDepth.field=refDepth.field, 
                 sampleNames.field=sampleNames.field)
    if (!all(strings <- sapply(args, isSingleString)))
        stop(paste0("'",names(args[!strings]), "'", collapse=","), 
             " must be character(1)")

    ## match fields
    matched.fields <- match(tolower(args), tolower(names(mcols(gr))), 0)
    found.fields <- matched.fields != 0L
    matched <- matched.fields[found.fields]

    ## error for 'ref' (required)
    if (!"ref.field" %in% names(args)[found.fields])
      stop("No 'ref' column could be identified.")

    ## extract fields, coerce to type
    lst <- as.list(c(rep(NA_character_, 2), rep(NA_integer_, 3), NA_character_))
    names(lst) <- names(args)
    type <- vapply(lst, class, character(1))
    lst[found.fields] <- Map(as, mcols(gr)[matched], type[found.fields])
    if (keep.extra.columns)
        extra <- mcols(gr)[, -matched, drop=FALSE]
    else
        extra <- DataFrame() 

    ## construct VRanges
    VRanges(seqnames=seqnames(gr), 
            seqinfo=seqinfo(gr),
            ranges=ranges(gr),
            ref=lst[["ref.field"]], 
            alt=lst[["alt.field"]], 
            totalDepth=lst[["totalDepth.field"]], 
            altDepth=lst[["altDepth.field"]], 
            refDepth=lst[["refDepth.field"]], 
            sampleNames=lst[["sampleNames.field"]], 
            extra)
}

setAs("GRanges", "VRanges",
    function(from) makeVRangesFromGRanges(from, keep.extra.columns=TRUE)
)
