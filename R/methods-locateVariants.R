### =========================================================================
### locateVariants methods 
### =========================================================================

### There are 7 defined variant regions :
### CodingVariants, IntronVariants, ThreeUTRVariants, FiveUTRVariants,
### IntergenicVariants, SpliceSiteVariants, AllVariants
### 
### Each variant region has methods for : 
###   query %in% Ranges, VCF, GRanges
###   subject %in% TranscriptDb, GRangesList 


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
        usub <- unlist(subject, use.names=FALSE)
        fo <- findOverlaps(query, usub, type="within")
        queryid <- queryHits(fo)
        if (length(fo) > 0) {
            txid <- rep(names(subject), elementLengths(subject))
            meta <- DataFrame(LOCATION=.location(length(queryid), "coding"),
                              QUERYID=queryid, 
                              TXID=as.integer(txid[subjectHits(fo)]), 
                              CDSID=values(usub)[["cds_id"]][subjectHits(fo)],
                              GENEID=rep(NA_character_, length(queryid)))
            .makeResult(query, meta)
        } else {
            .makeResult(query, NULL)
        } 
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
    genedf <- 
      data.frame(geneid=rep(names(cache[["txbygene"]]),
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
        .makeResult(query, .makeMeta(query, subject, "intron"))
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
        values(res)[["geneID"]] <-
            genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "ThreeUTRVariants"),
    function(query, subject, region, ...)
    {
        .makeResult(query, .makeMeta(query, subject, "threeUTR"))
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
        values(res)[["geneID"]] <-
            genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "FiveUTRVariants"),
    function(query, subject, region, ...)
{
    .makeResult(query, .makeMeta(query, subject, "fiveUTR"))
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
    .makeResult(query, .intergenic(query, subject))
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
        values(res)[["geneID"]] <-
            genedf$geneid[match(values(res)[["TXID"]], genedf$txid)]
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", 
          "SpliceSiteVariants"),
    function(query, subject, region, ...)
    {
        .makeResult(query, .spliceSites(query, subject))
    }
)

### -------------------------------------------------------------------------
### region = AllVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "AllVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()))
    {
        coding <- locateVariants(query, subject, 
                                 CodingVariants(), cache=cache)
        intron <- locateVariants(query, subject, 
                                 IntronVariants(), cache=cache)
        splice <- locateVariants(query, subject, 
                                 SpliceSiteVariants(), cache=cache)
        fiveUTR <- locateVariants(query, subject, 
                                  FiveUTRVariants(), cache=cache)
        threeUTR <- locateVariants(query, subject, 
                                   ThreeUTRVariants(), cache=cache)
        intergenic <- locateVariants(query, subject, 
                                     IntergenicVariants(), cache=cache)

        base <- c(coding, intron, fiveUTR, threeUTR, splice)
        precedesID <- followsID <- rep(NA_character_, length(base))
        values(base) <- append(values(base), DataFrame(precedesID, followsID))
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
    int_start <- GRanges(seqnames(usub),
                         IRanges(start(usub), start(usub) + 1),
                         strand=strand(usub))
    int_end <- GRanges(seqnames(usub),
                       IRanges(end(usub) - 1, end(usub)),
                       strand=strand(usub))
    fo_start <- findOverlaps(query, int_start, type="any")
    fo_end <- findOverlaps(query, int_end, type="any")
    ## FIXME : no method for 'c' for Hits class
    fo <- data.frame(rbind(as.matrix(fo_start), as.matrix(fo_end)))
    fo <- fo[!duplicated(fo),]
    queryid <- fo$queryHits
    if (length(fo) > 0) {
        txid <- rep(names(subject), elementLengths(subject))
        meta <- DataFrame(LOCATION=.location(length(queryid), "spliceSite"),
                          QUERYID=queryid, 
                          TXID=as.integer(txid[fo$subjectHits]), 
                          CDSID=rep(NA_integer_, length(queryid)),
                          GENEID=rep(NA_character_, length(queryid)))
    } else {
        NULL
    }
}

.intergenic <- function(query, subject, ...)
{
    usub <- range(subject)
    co <- countOverlaps(query, subject, type="any")
    intergenic <- co == 0
    if (all(intergenic == FALSE)) {
        DataFrame(LOCATION=.location(), QUERYID=integer(),
                  TXID=integer(), CDSID=integer(), GENEID=character(),
                  PRECEDEID=character(), FOLLOWID=character())
    } else {
        query <- query[intergenic]

        ## gene ranges
        rnglst <- subject
        rng <- unlist(rnglst, use.names=FALSE)
        genes <- rep(names(rnglst), elementLengths(rnglst))
        ## query precedes subject; get index for following gene 
        pidx <- precede(query, rng)
        ## query follows subject; get index for preceding gene 
        fidx <- follow(query, rng)

        location <- .location(length(query))
        location[as.vector(seqnames(query)) %in% seqlevels(rng)] <- "intergenic"

        DataFrame(LOCATION=location, QUERYID=which(intergenic),
                  TXID=NA_integer_, CDSID=NA_integer_,
                  GENEID=NA_character_,
                  PRECEDEID=genes[pidx], FOLLOWID=genes[fidx])
    }
}

.makeMeta <- function(query, subject, vtype, ...)
{
    ## elementMetadata are not propagated
    ## columns must match for combining in AllVariants
    usub <- unlist(subject, use.names=FALSE)
    fo <- findOverlaps(query, usub, type="within")
    queryid <- queryHits(fo)
    if (length(fo) > 0) {
        txid <- rep(names(subject), elementLengths(subject))
        DataFrame(LOCATION=.location(length(queryid), vtype),
                  QUERYID=queryid, 
                  TXID=as.integer(txid[subjectHits(fo)]), 
                  CDSID=rep(NA_integer_, length(queryid)),
                  GENEID=rep(NA_character_, length(queryid)))
    } else {
        NULL
    }
}

.makeResult <- function(query, meta,  ...)
{
    if (!is.null(meta)) {
        if (0L == nrow(meta)) {
            res <- GRanges()
            values(res) <- meta
        } else {
            res <- query[meta$QUERYID]
            names(res) <- names(query)[meta$QUERYID]
            values(res) <- meta
            res
       }
    } else {
        res <- GRanges()
        values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(),
                                 TXID=integer(), CDSID=integer(),
                                 GENEID=character())
    }
    res
}

