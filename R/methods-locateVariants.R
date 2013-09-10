### =========================================================================
### locateVariants methods 
### =========================================================================

### The 8 defined variant regions :
### CodingVariants, IntronVariants, ThreeUTRVariants, FiveUTRVariants,
### IntergenicVariants, SpliceSiteVariants, PromoterVariants, AllVariants
### 
### Each region has the following methods : 
### query %in% Ranges, VCF, GRanges
### subject %in% TranscriptDb, GRangesList 


### -------------------------------------------------------------------------
### methods applicable to all variant regions 
###

setMethod("locateVariants", c("Ranges", "TranscriptDb", "VariantType"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
{
    callGeneric(as(query, "GRanges"), subject, region, ..., cache=cache,
        ignore.strand=ignore.strand, asHits=asHits)
})

setMethod("locateVariants", c("Ranges", "GRangesList", "VariantType"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
{
    callGeneric(as(query, "GRanges"), subject, region, ..., cache=cache,
        ignore.strand=ignore.strand, asHits=asHits)
})

setMethod("locateVariants", c("VCF", "TranscriptDb", "VariantType"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
{
    callGeneric(rowData(query), subject, region, ..., cache=cache,
        ignore.strand=ignore.strand, asHits=asHits)
})

setMethod("locateVariants", c("VCF", "GRangesList", "VariantType"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
{
    callGeneric(rowData(query), subject, region, ..., cache=cache,
        ignore.strand=ignore.strand, asHits=asHits)
})

### -------------------------------------------------------------------------
### region = CodingVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "CodingVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
    {
        if (!exists("mask", cache, inherits=FALSE)) {
            origSeqlevels <- seqlevels(subject)
            if (.setSubjectSeq(query, subject)) 
                on.exit(seqlevels(subject) <- origSeqlevels)
        } 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1
        if (!exists("cdsbytx", cache, inherits=FALSE))
            cache[["cdsbytx"]] <- cdsBy(subject)
        res <- callGeneric(query, cache[["cdsbytx"]], region, ..., 
            ignore.strand=ignore.strand, asHits=asHits)
        if (is(res, "GenomicRanges") & length(res) > 0L) {
            res$GENEID <- select(subject, res$TXID, "GENEID", "TXID")$GENEID 
        }
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "CodingVariants"),
    function(query, subject, region, ..., ignore.strand=FALSE, asHits=FALSE)
    {
        .makeResult(query, subject, "coding", ignore.strand=ignore.strand,
            asHits=asHits)
    }
)

### -------------------------------------------------------------------------
### region = IntronVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "IntronVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
    {
        if (!exists("mask", cache, inherits=FALSE)) {
            origSeqlevels <- seqlevels(subject)
            if (.setSubjectSeq(query, subject)) 
                on.exit(seqlevels(subject) <- origSeqlevels)
        } 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1
        if (!exists("intbytx", cache, inherits=FALSE))
            cache[["intbytx"]] <- intronsByTranscript(subject)
        res <- callGeneric(query, cache[["intbytx"]], region, ...,
            ignore.strand=ignore.strand, asHits=asHits)
        if (is(res, "GenomicRanges") & length(res) > 0L) {
            res$GENEID <- select(subject, res$TXID, "GENEID", "TXID")$GENEID 
        }
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "IntronVariants"),
    function(query, subject, region, ..., ignore.strand=FALSE, asHits=FALSE)
    {
        .makeResult(query, subject, "intron", ignore.strand=ignore.strand,
            asHits=asHits)
    }
)

### -------------------------------------------------------------------------
### region = ThreeUTRVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "ThreeUTRVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
    {
        if (!exists("mask", cache, inherits=FALSE)) {
            origSeqlevels <- seqlevels(subject)
            if (.setSubjectSeq(query, subject)) 
                on.exit(seqlevels(subject) <- origSeqlevels)
        } 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0)) 
            start(query)[insertion] <- start(query)[insertion] - 1 

        if (!exists("threeUTRbytx", cache, inherits=FALSE))
            cache[["threeUTRbytx"]] <- threeUTRsByTranscript(subject)
        res <- callGeneric(query, cache[["threeUTRbytx"]], region, ...,
            ignore.strand=ignore.strand, asHits=asHits)
        if (is(res, "GenomicRanges") & length(res) > 0L) {
            res$GENEID <- select(subject, res$TXID, "GENEID", "TXID")$GENEID 
        }
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "ThreeUTRVariants"),
    function(query, subject, region, ..., ignore.strand=FALSE, asHits=FALSE)
    {
        .makeResult(query, subject, "threeUTR", ignore.strand=ignore.strand,
                    asHits=asHits)
    }
)

### -------------------------------------------------------------------------
### region = FiveUTRVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "FiveUTRVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
    {
        if (!exists("mask", cache, inherits=FALSE)) {
            origSeqlevels <- seqlevels(subject)
            if (.setSubjectSeq(query, subject)) 
                on.exit(seqlevels(subject) <- origSeqlevels)
        } 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1

        if (!exists("fiveUTRbytx", cache, inherits=FALSE))
            cache[["fiveUTRbytx"]] <- fiveUTRsByTranscript(subject)
        res <- callGeneric(query, cache[["fiveUTRbytx"]], region, ...,
            ignore.strand=ignore.strand, asHits=asHits)
        if (is(res, "GenomicRanges") & length(res) > 0L) {
            res$GENEID <- select(subject, res$TXID, "GENEID", "TXID")$GENEID 
        }
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "FiveUTRVariants"),
    function(query, subject, region, ..., ignore.strand=FALSE, asHits=FALSE)
    {
        .makeResult(query, subject, "fiveUTR", ignore.strand=ignore.strand,
            asHits=asHits)
    }
)

### -------------------------------------------------------------------------
### region = IntergenicVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", 
          "IntergenicVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE)
    {
        if (!exists("mask", cache, inherits=FALSE)) {
            origSeqlevels <- seqlevels(subject)
            if (.setSubjectSeq(query, subject)) 
                on.exit(seqlevels(subject) <- origSeqlevels)
        } 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1

        if (!exists("txbygene", cache, inherits=FALSE))
            cache[["txbygene"]] <- transcriptsBy(subject, "gene")

        callGeneric(query, cache[["txbygene"]], region, ..., 
            ignore.strand=ignore.strand)
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "IntergenicVariants"),
    function(query, subject, region, ..., ignore.strand=FALSE)
    {
        .intergenic(query, subject, region, ignore.strand=ignore.strand)
    }
)

### -------------------------------------------------------------------------
## region = SpliceSiteVariants 
##

setMethod("locateVariants", c("GRanges", "TranscriptDb", "SpliceSiteVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
    {
        if (!exists("mask", cache, inherits=FALSE)) {
            origSeqlevels <- seqlevels(subject)
            if (.setSubjectSeq(query, subject)) 
                on.exit(seqlevels(subject) <- origSeqlevels)
        } 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1

        if (!exists("intbytx", cache, inherits=FALSE))
            cache[["intbytx"]] <- intronsByTranscript(subject)
        res <- callGeneric(query, cache[["intbytx"]], region, ...,
            ignore.strand=ignore.strand, asHits=asHits)
        if (is(res, "GenomicRanges") & length(res) > 0L) {
            res$GENEID <- select(subject, res$TXID, "GENEID", "TXID")$GENEID 
        }
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "SpliceSiteVariants"),
    function(query, subject, region, ..., ignore.strand=FALSE, asHits=FALSE)
    {
        .spliceSites(query, subject, ignore.strand=ignore.strand, asHits=asHits)
    }
)

### -------------------------------------------------------------------------
### region = PromoterVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "PromoterVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE, asHits=FALSE)
    { 
        if (!exists("mask", cache, inherits=FALSE)) {
            origSeqlevels <- seqlevels(subject)
            if (.setSubjectSeq(query, subject)) 
                on.exit(seqlevels(subject) <- origSeqlevels)
        } 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())

        ## for width(ranges) == 0 : de-increment start to equal end value 
        if (any(insertion <- width(query) == 0)) 
            start(query)[insertion] <- start(query)[insertion] - 1 

        if (!exists("tx", cache, inherits=FALSE)) {
            tx <- transcripts(subject) 
            cache[["tx"]] <- splitAsList(tx, seq_len(length(tx))) 
            names(cache[["tx"]]) <- tx$tx_id
        }
        res <- callGeneric(query, cache[["tx"]], region, ...,
            ignore.strand=ignore.strand, asHits=asHits)
        if (is(res, "GenomicRanges") & length(res) > 0L) {
            if (!is.null(res$TXID))
                res$GENEID <- select(subject, res$TXID, "GENEID", "TXID")$GENEID
        }
        res
    }
)

setMethod("locateVariants", c("GRanges", "GRangesList", "PromoterVariants"),
    function(query, subject, region, ..., ignore.strand=FALSE, asHits=FALSE)
    {

        if (is.null(txid <- names(subject)))
            txid <- NA_integer_
        usub <- unlist(subject, use.names=FALSE)
        pm <- promoters(usub, upstream(region), downstream(region))
        fo <- findOverlaps(query, pm, type="within", 
            ignore.strand=ignore.strand)

        if (asHits)
            return(.consolidateHits(fo, length(query), length(subject),
                elementLengths(subject)))

        if (length(fo) > 0) {
            queryid <- queryHits(fo)
            GRanges(seqnames=seqnames(query)[queryid],
                    ranges=IRanges(ranges(query)[queryid]),
                    strand=strand(query)[queryid],
                    LOCATION=.location(length(queryid), "promoter"), 
                    QUERYID=queryid,
                    TXID=as.integer(txid[subjectHits(fo)]),
                    CDSID=NA_integer_,
                    GENEID=NA_character_,
                    PRECEDEID=CharacterList(character(0)),
                    FOLLOWID=CharacterList(character(0)))
        } else {
            res <- GRanges()
            values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(),
                                     TXID=integer(), CDSID=integer(),
                                     GENEID=character(),
                                     PRECEDEID=CharacterList(),
                                     FOLLOWID=CharacterList())
            res
        }
    }
)

### -------------------------------------------------------------------------
### region = AllVariants 
###

setMethod("locateVariants", c("GRanges", "TranscriptDb", "AllVariants"),
    function(query, subject, region, ..., cache=new.env(parent=emptyenv()),
             ignore.strand=FALSE)
    {
        origSeqlevels <- seqlevels(subject)
        if (.setSubjectSeq(query, subject)) 
            on.exit(seqlevels(subject) <- origSeqlevels) 
        if (!any(seqlevels(query) %in% seqlevels(subject)))
            return(.returnEmpty())
        cache[["mask"]] <- TRUE

        coding <- 
            locateVariants(query, subject, CodingVariants(), cache=cache,
                           ignore.strand=ignore.strand)
        intron <- 
            locateVariants(query, subject, IntronVariants(), cache=cache,
                           ignore.strand=ignore.strand)
        splice <- 
            locateVariants(query, subject, SpliceSiteVariants(), 
                           cache=cache, ignore.strand=ignore.strand)
        promoter <- 
            locateVariants(query, subject,
                           PromoterVariants(upstream(promoter(region)),
                                            downstream(promoter(region))), 
                           cache=cache, ignore.strand=ignore.strand)

        ## Consolidate calls for UTR data:
        if (!exists("fiveUTRbytx", cache, inherits=FALSE)) {
            splicings <- 
                GenomicFeatures:::.getSplicingsForTranscriptsWithCDSs(subject)
            cache[["fiveUTRbytx"]] <- 
                GenomicFeatures:::.make5UTRsByTranscript(subject, splicings)
        }
        if (!exists("threeUTRbytx", cache, inherits=FALSE))
            cache[["threeUTRbytx"]] <- 
                GenomicFeatures:::.make3UTRsByTranscript(subject, splicings)

        fiveUTR <- 
            locateVariants(query, subject, FiveUTRVariants(), 
                           cache=cache, ignore.strand=ignore.strand)
        threeUTR <- 
            locateVariants(query, subject, ThreeUTRVariants(), 
                           cache=cache, ignore.strand=ignore.strand)
        intergenic <- 
            locateVariants(query, subject, 
                           IntergenicVariants(upstream(intergenic(region)),
                                              downstream(intergenic(region))), 
                           cache=cache, ignore.strand=ignore.strand)

        ans <- c(coding, intron, fiveUTR, threeUTR, splice, promoter, intergenic)
        meta <- values(ans)
        ans[order(meta$QUERYID, meta$TXID, meta$CDSID, meta$GENEID), ]
    }
)

### -------------------------------------------------------------------------
### helpers 
###

.returnEmpty <- function()
{
    warning("none of seqlevels(query) match seqlevels(subject)")
    res <- GRanges()
    values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(), 
                             TXID=integer(), CDSID=integer(), 
                             GENEID=character(), PRECEDEID=CharacterList(), 
                             FOLLOWID=CharacterList())
    res
}

.location <-
    function(length=0, value=NA)
{
    levels <- c("spliceSite", "intron", "fiveUTR", "threeUTR",
        "coding", "intergenic", "promoter")
    factor(rep(value, length), levels=levels)
}

.spliceSites <- function(query, subject, ignore.strand, asHits, ...)
{
    ## Overlap any portion of first 2 and last 2 nucleotides of introns.
    usub <- unlist(subject, use.names=FALSE)
    int_start <- GRanges(seqnames(usub), IRanges(start(usub), start(usub) + 1),
        strand=strand(usub))
    int_end <- GRanges(seqnames(usub), IRanges(end(usub) - 1, end(usub)),
        strand=strand(usub))
    fo_start <- findOverlaps(query, int_start, type="any",
        ignore.strand=ignore.strand)
    fo_end <- findOverlaps(query, int_end, type="any",
        ignore.strand=ignore.strand)
    fo <- union(fo_start, fo_end)

    if (asHits)
        return(.consolidateHits(fo, length(query), length(subject),
            elementLengths(subject))) 
    if (length(fo) > 0) {
        txid <- rep(names(subject), elementLengths(subject))
        queryid <- queryHits(fo)
        GRanges(seqnames=seqnames(query)[queryid],
                ranges=IRanges(ranges(query)[queryid]),
                strand=strand(query)[queryid],
                LOCATION=.location(length(queryid), "spliceSite"),
                QUERYID=queryid,
                TXID=as.integer(txid[subjectHits(fo)]),
                CDSID=NA_integer_,
                GENEID=NA_character_,
                PRECEDEID=CharacterList(character(0)),
                FOLLOWID=CharacterList(character(0)))
    } else {
        res <- GRanges()
        values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(),
                                 TXID=integer(), CDSID=integer(), 
                                 GENEID=character(), PRECEDEID=CharacterList(),
                                 FOLLOWID=CharacterList())
        res
    }
}

.consolidateHits <- function(hits, qlen, slen, elen, ...)
{
    txid <- rep(seq_len(slen), elen)
    dat <- cbind(queryHits(hits), as.integer(txid[subjectHits(hits)]))
    unq <- dat[!duplicated(dat),]
    return(new("Hits", 
           queryHits=unq[,1], 
           subjectHits=unq[,2],
           queryLength=qlen, 
           subjectLength=slen)) 
}

.intergenic <- function(query, subject, region, ignore.strand, ...)
{
    co <- countOverlaps(query, subject, type="any", 
                        ignore.strand=ignore.strand)
    intergenic <- co == 0
    if (all(!intergenic | (length(subject) == 0L))) {
        res <- GRanges()
        values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(), 
                                 TXID=integer(), CDSID=integer(), 
                                 GENEID=character(), PRECEDEID=CharacterList(), 
                                 FOLLOWID=CharacterList())
        res
    } else {
        q_range <- query[intergenic]
        s_range <- unlist(subject, use.names=FALSE)
        s_genes <- rep(names(subject), elementLengths(subject))

        ## ID all genes that fall in upstream / downstream range.
        ## upstream == follow:
        f_range <- .shiftRangeUpDown(q_range, upstream(region), TRUE)
        f_fo <- findOverlaps(f_range, s_range)
        f_genes <- unname(splitAsList(s_genes[subjectHits(f_fo)],
                                factor(queryHits(f_fo),
                                       levels=seq_len(queryLength(f_fo)))))
        f_genes <- unique(f_genes)
 
        ## downstream == precede:
        p_range <- .shiftRangeUpDown(q_range, downstream(region), FALSE)
        p_fo <- findOverlaps(p_range, s_range)
        p_genes <- unname(splitAsList(s_genes[subjectHits(p_fo)],
                                factor(queryHits(p_fo),
                                       levels=seq_len(queryLength(p_fo)))))
        p_genes <- unique(p_genes)

        values(q_range) <- 
            DataFrame(LOCATION=.location(length(q_range), "intergenic"),
                      QUERYID=which(intergenic), TXID=NA_integer_, 
                      CDSID=NA_integer_, GENEID=NA_character_, 
                      PRECEDEID=p_genes, FOLLOWID=f_genes)
        q_range 
    }
}

## Behaves like flank(); ranges do not include start/end values.
.shiftRangeUpDown <- function(x, distance, upstream=TRUE)
{
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_minus <- which(strand(x) == "-")
    if (upstream) {
        ## '+' strand
        on_plus_end <- start(x)[on_plus] - 1L 
        on_plus_start <- start(x)[on_plus] - distance 
        ## '-' strand
        on_minus_start <- end(x)[on_minus] + 1L 
        on_minus_end <- end(x)[on_minus] + distance 
    } else {
        ## '+' strand
        on_plus_start <- end(x)[on_plus] + 1L 
        on_plus_end <- end(x)[on_plus] + distance
        ## '-' strand 
        on_minus_start <- start(x)[on_minus] - distance
        on_minus_end <- start(x)[on_minus] - 1L 
    }
    start(x)[on_plus] <- on_plus_start
    end(x)[on_plus] <- on_plus_end
    start(x)[on_minus] <- on_minus_start
    end(x)[on_minus] <- on_minus_end
    x
}

.makeResult <- function(query, subject, vtype, ignore.strand, asHits)
{
    if (asHits) {
        return(findOverlaps(query, subject, type="within", 
            ignore.strand=ignore.strand))
    } else {
        usub <- unlist(subject, use.names=FALSE)
        fo <- findOverlaps(query, usub, type="within", 
            ignore.strand=ignore.strand)
    }
    if (length(fo) > 0) {
        queryid <- queryHits(fo)
        txid <- NA_integer_
        cdsid <- NA_integer_
        if (!is.null(tx <- rep(names(subject), elementLengths(subject))))
            txid <- tx 
        if (vtype == "coding")
            if(!is.null(cds <- values(usub)[["cds_id"]][subjectHits(fo)]))
                cdsid <- cds

        GRanges(seqnames=seqnames(query)[queryid],
                ranges=IRanges(ranges(query)[queryid]),
                strand=strand(query)[queryid],
                LOCATION=.location(length(queryid), vtype),
                QUERYID=queryid,
                TXID=as.integer(txid[subjectHits(fo)]),
                CDSID=cdsid,
                GENEID=NA_character_,
                PRECEDEID=CharacterList(character(0)),
                FOLLOWID=CharacterList(character(0)))
    } else {
        res <- GRanges()
        values(res) <- DataFrame(LOCATION=.location(), QUERYID=integer(),
                                 TXID=integer(), CDSID=integer(), 
                                 GENEID=character(), PRECEDEID=CharacterList(),
                                 FOLLOWID=CharacterList())
        res
    }
}
