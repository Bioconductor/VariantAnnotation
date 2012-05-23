### =========================================================================
### locateVariants methods 
### =========================================================================

### The 7 defined variant regions :
### CodingVariants, IntronVariants, ThreeUTRVariants, FiveUTRVariants,
### IntergenicVariants, SpliceSiteVariants, AllVariants
### 
### Each variant region has the following methods : 
### query %in% Ranges, VCF, GRanges
### subject %in% TranscriptDb, GRangesList 


### -------------------------------------------------------------------------
### methods applicable to all variant regions 
###

setMethod("locateVariants", c("Ranges", "TranscriptDb", "ANY"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
{
    callGeneric(as(query, "GRanges"), subject, region, ..., cache=cache)
})
setMethod("locateVariants", c("Ranges", "GRangesList", "ANY"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
{
    callGeneric(as(query, "GRanges"), subject, region, ..., cache=cache)
})

setMethod("locateVariants", c("VCF", "TranscriptDb", "ANY"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
{
    callGeneric(rowData(query), subject, region, 
                ..., cache=cache)
})
setMethod("locateVariants", c("VCF", "GRangesList", "ANY"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
{
    callGeneric(rowData(query), subject, region, 
                ..., cache=cache)
})

### -------------------------------------------------------------------------
## region = CodingVariants 
##

setMethod("locateVariants", c("GRanges", "TranscriptDb", "CodingVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
{
    queryseq <- seqlevels(query)
    subseq <- seqlevels(subject)
    if (!any(queryseq %in% subseq))
        warning("none of seqlevels(query) match seqlevels(subject)")

    ## mask chromosomes not in query
    masks <- isActiveSeq(subject)
    on.exit(isActiveSeq(subject) <- masks)
    .setActiveSubjectSeq(query, subject)

    ## for width(ranges) == 0 : de-increment start to equal end value 
    if (any(insertion <- width(query) == 0))
        start(query)[insertion] <- start(query)[insertion] - 1

    if (!exists("cdsbytx", cache, inherits=FALSE))
        cache[["cdsbytx"]] <- cdsBy(subject)
    if (!exists("txbygene", cache, inherits=FALSE))
        cache[["txbygene"]] <- transcriptsBy(subject, "gene")

    res <- callGeneric(query, cache[["cdsbytx"]], region, ...)
    genedf <- data.frame(geneid=rep(names(cache[["txbygene"]]), 
                             elementLengths(cache[["txbygene"]])),
                         txid=values(unlist(cache[["txbygene"]], 
                             use.names=FALSE))[["tx_id"]], 
                             stringsAsFactors=FALSE)
    values(res)[["GENEID"]] <-
        genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
    res
})

setMethod("locateVariants", c("GRanges", "GRangesList", "CodingVariants"),
    function(query, subject, region, ...)
    {
        .makeResult(query, subject, "coding")
    }
)

### -------------------------------------------------------------------------
### region = IntronVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "IntronVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
{
    queryseq <- seqlevels(query)
    subseq <- seqlevels(subject)
    if (!any(queryseq %in% subseq))
        warning("none of seqlevels(query) match seqlevels(subject)")

    ## mask chromosomes not in query
    masks <- isActiveSeq(subject)
    on.exit(isActiveSeq(subject) <- masks)
    .setActiveSubjectSeq(query, subject)

    ## for width(ranges) == 0 : de-increment start to equal end value 
    if (any(insertion <- width(query) == 0))
        start(query)[insertion] <- start(query)[insertion] - 1

    if (!exists("intbytx", cache, inherits=FALSE))
        cache[["intbytx"]] <- intronsByTranscript(subject)
    if (!exists("txbygene", cache, inherits=FALSE))
        cache[["txbygene"]] <- transcriptsBy(subject, "gene")

    res <- callGeneric(query, cache[["intbytx"]], region, ...)
    genedf <- data.frame(geneid=rep(names(cache[["txbygene"]]),
                             elementLengths(cache[["txbygene"]])),
                         txid=values(unlist(cache[["txbygene"]], 
                             use.names=FALSE))[["tx_id"]],
                             stringsAsFactors=FALSE)
    values(res)[["GENEID"]] <- 
        genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
    res
})


setMethod("locateVariants", c("GRanges", "GRangesList", "IntronVariants"),
    function(query, subject, region, ...)
    {
        .makeResult(query, subject, "intron")
    }
)

### -------------------------------------------------------------------------
### region = ThreeUTRVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "ThreeUTRVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
    {
        queryseq <- seqlevels(query)
        subseq <- seqlevels(subject)
        if (!any(queryseq %in% subseq))
            warning("none of seqlevels(query) match seqlevels(subject)")

        ## mask chromosomes not in query
        masks <- isActiveSeq(subject)
        on.exit(isActiveSeq(subject) <- masks)
        .setActiveSubjectSeq(query, subject)

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0)) 
            start(query)[insertion] <- start(query)[insertion] - 1 

        if (!exists("threeUTRbytx", cache, inherits=FALSE))
            cache[["threeUTRbytx"]] <- threeUTRsByTranscript(subject)
        if (!exists("txbygene", cache, inherits=FALSE))
            cache[["txbygene"]] <- transcriptsBy(subject, "gene")

        res <- callGeneric(query, cache[["threeUTRbytx"]], region, ...)
        genedf <- data.frame(geneid = rep(names(cache[["txbygene"]]),
                                 elementLengths(cache[["txbygene"]])),
                             txid = values(unlist(cache[["txbygene"]],
                                 use.names=FALSE))[["tx_id"]], 
                                 stringsAsFactors=FALSE)
        values(res)[["GENEID"]] <-
            genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "ThreeUTRVariants"),
    function(query, subject, region, ...)
    {
        .makeResult(query, subject, "threeUTR")
    }
)

### -------------------------------------------------------------------------
### region = FiveUTRVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "FiveUTRVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
    {
        queryseq <- seqlevels(query)
        subseq <- seqlevels(subject)
        if (!any(queryseq %in% subseq))
            warning("none of seqlevels(query) match seqlevels(subject)")

        ## mask chromosomes not in query
        masks <- isActiveSeq(subject)
        on.exit(isActiveSeq(subject) <- masks)
        .setActiveSubjectSeq(query, subject)

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1

        if (!exists("fiveUTRbytx", cache, inherits=FALSE))
            cache[["fiveUTRbytx"]] <- fiveUTRsByTranscript(subject)
        if (!exists("txbygene", cache, inherits=FALSE))
            cache[["txbygene"]] <- transcriptsBy(subject, "gene")

        res <- callGeneric(query, cache[["fiveUTRbytx"]], region, ...)
        genedf <- data.frame(geneid = rep(names(cache[["txbygene"]]),
                                 elementLengths(cache[["txbygene"]])),
                             txid = values(unlist(cache[["txbygene"]],
                                 use.names=FALSE))[["tx_id"]],
                                 stringsAsFactors=FALSE)
        values(res)[["GENEID"]] <-
            genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "FiveUTRVariants"),
    function(query, subject, region, ...)
{
    .makeResult(query, subject, "fiveUTR")
})

### -------------------------------------------------------------------------
### region = IntergenicVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", 
          "IntergenicVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
    {
        queryseq <- seqlevels(query)
        subseq <- seqlevels(subject)
        if (!any(queryseq %in% subseq))
            warning("none of seqlevels(query) match seqlevels(subject)")

        ## mask chromosomes not in query
        masks <- isActiveSeq(subject)
        on.exit(isActiveSeq(subject) <- masks)
        .setActiveSubjectSeq(query, subject)

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1

        if (!exists("txbygene", cache, inherits=FALSE))
            cache[["txbygene"]] <- transcriptsBy(subject, "gene")

        callGeneric(query, cache[["txbygene"]], region, ...)
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "IntergenicVariants"),
    function(query, subject, region, ...)
{
    .intergenic(query, subject)
})

### -------------------------------------------------------------------------
## region = SpliceSiteVariants 
##

setMethod("locateVariants", c("GRanges", "TranscriptDb",
          "SpliceSiteVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
    {
        queryseq <- seqlevels(query)
        subseq <- seqlevels(subject)
        if (!any(queryseq %in% subseq))
            warning("none of seqlevels(query) match seqlevels(subject)")

        ## mask chromosomes not in query
        masks <- isActiveSeq(subject)
        on.exit(isActiveSeq(subject) <- masks)
        .setActiveSubjectSeq(query, subject)

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1

        if (!exists("intbytx", cache, inherits=FALSE))
            cache[["intbytx"]] <- intronsByTranscript(subject)
        if (!exists("txbygene", cache, inherits=FALSE))
            cache[["txbygene"]] <- transcriptsBy(subject, "gene")

        res <- callGeneric(query, cache[["intbytx"]], region, ...)
        genedf <- data.frame(geneid=rep(names(cache[["txbygene"]]),
                                 elementLengths(cache[["txbygene"]])),
                             txid=values(unlist(cache[["txbygene"]],
                                 use.names=FALSE))[["tx_id"]],
                                 stringsAsFactors=FALSE)
        values(res)[["GENEID"]] <-
            genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", 
          "SpliceSiteVariants"),
    function(query, subject, region, ...)
    {
        .spliceSites(query, subject)
    }
)

### -------------------------------------------------------------------------
### region = AllVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "AllVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
    {
        coding <- locateVariants(query, subject, CodingVariants(), cache=cache)
        intron <- locateVariants(query, subject, IntronVariants(), cache=cache)
        splice <- locateVariants(query, subject, SpliceSiteVariants(), 
                                 cache=cache)

        ## Consolidate calls for UTR data
        if (!exists("fiveUTRbytx", cache, inherits=FALSE)) {
            splicings <- 
                GenomicFeatures:::.getSplicingsForTranscriptsWithCDSs(subject)
            cache[["fiveUTRbytx"]] <- 
                GenomicFeatures:::.make5UTRsByTranscript(txdb, splicings)
        }
        if (!exists("threeUTRbytx", cache, inherits=FALSE))
            cache[["threeUTRbytx"]] <- 
                GenomicFeatures:::.make5UTRsByTranscript(txdb, splicings)

        fiveUTR <- locateVariants(query, subject, FiveUTRVariants(), 
                                  cache=cache)
        threeUTR <- locateVariants(query, subject, ThreeUTRVariants(), 
                                   cache=cache)
        intergenic <- locateVariants(query, subject, IntergenicVariants(), 
                                     cache=cache)

        base <- c(coding, intron, fiveUTR, threeUTR, splice)
        PRECEDEID <- FOLLOWID <- rep(NA_character_, length(base))
        values(base) <- append(values(base), DataFrame(PRECEDEID, FOLLOWID))
        ans <- c(base, intergenic)
        meta <- values(ans)
        ans[order(meta$QUERYID, meta$TXID, meta$CDSID, meta$GENEID), ]
    }
)

### -------------------------------------------------------------------------
### helpers 
###

.location <-
    function(length=0, value=NA)
{
    levels <- c("spliceSite", "intron", "fiveUTR", "threeUTR",
        "coding", "intergenic")
    factor(rep(value, length), levels=levels)
}

.spliceSites <- function(query, subject, ...)
{
    ## overlap any portion of first 2 and last 2 nucleotides of introns
    usub <- unlist(subject, use.names=FALSE)
    int_start <- GRanges(seqnames(usub), IRanges(start(usub), start(usub) + 1),
        strand=strand(usub))
    int_end <- GRanges(seqnames(usub), IRanges(end(usub) - 1, end(usub)),
        strand=strand(usub))
    fo_start <- findOverlaps(query, int_start, type="any")
    fo_end <- findOverlaps(query, int_end, type="any")
    fo <- union(fo_start, fo_end)

    if (length(fo) > 0) {
        queryid <- queryHits(fo)
        txid <- rep(names(subject), elementLengths(subject))

        GRanges(seqnames=seqnames(query)[queryid],
                ranges=IRanges(ranges(query)[queryid]),
                strand=strand(query)[queryid],
                LOCATION=.location(length(queryid), "spliceSite"),
                QUERYID=queryid,
                TXID=as.integer(txid[subjectHits(fo)]),
                CDSID=NA_integer_,
                GENEID=NA_character_)
    } else {
        res <- GRanges()
        values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(),
                                 TXID=integer(), CDSID=integer(), 
                                 GENEID=character())
        res
    }
}

.intergenic <- function(query, subject, ...)
{
    co <- countOverlaps(query, subject, type="any")
    intergenic <- co == 0
    if (all(!intergenic)) {
        res <- GRanges()
        values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(), 
                                 TXID=integer(), CDSID=integer(), 
                                 GENEID=character(), PRECEDEID=character(), 
                                 FOLLOWID=character())
        res
    } else {
        res <- query[intergenic]

        ## gene ranges
        rng <- unlist(subject, use.names=FALSE)
        genes <- rep(names(subject), elementLengths(subject))
        ## query precedes subject; get index for following gene 
        pidx <- precede(res, rng)
        ## query follows subject; get index for preceding gene 
        fidx <- follow(res, rng)

        values(res) <- DataFrame(LOCATION=.location(length(res), "intergenic"),
                                 QUERYID=which(intergenic), TXID=NA_integer_, 
                                 CDSID=NA_integer_, GENEID=NA_character_, 
                                 PRECEDEID=genes[pidx], FOLLOWID=genes[fidx])
        res
    }
}

.makeResult <- function(query, subject, vtype, ...)
{
    usub <- unlist(subject, use.names=FALSE)
    fo <- findOverlaps(query, usub, type="within")
    if (length(fo) > 0) {
        queryid <- queryHits(fo)
        txid <- rep(names(subject), elementLengths(subject))
        if (vtype == "coding")
            cdsid <- values(usub)[["cds_id"]][subjectHits(fo)]
        else
            cdsid <- NA_integer_

        GRanges(seqnames=seqnames(query)[queryid],
                ranges=IRanges(ranges(query)[queryid]),
                strand=strand(query)[queryid],
                LOCATION=.location(length(queryid), vtype),
                QUERYID=queryid,
                TXID=as.integer(txid[subjectHits(fo)]),
                CDSID=cdsid,
                GENEID=NA_character_)
    } else {
        res <- GRanges()
        values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(),
                                 TXID=integer(), CDSID=integer(), 
                                 GENEID=character())
        res
    }
}

