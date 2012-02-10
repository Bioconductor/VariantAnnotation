### =========================================================================
### locateVariants methods 
### =========================================================================

setMethod("locateVariants",  c(query="Ranges", subject="TranscriptDb"), 
    function(query, subject, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject, ...) 
    }
)

setMethod("locateVariants", c(query="VCF", subject="TranscriptDb"), 
    function(query, subject, ...)
    {
        callGeneric(query=rowData(query), subject, ...) 
    }
)

.location <-
    function(length=0, value=NA)
{
    levels <- c("transcript_region", "intron", "5'UTR", "3'UTR",
        "coding", "intergenic")
    factor(rep(value, length), levels=levels)
}

setMethod("locateVariants", c("GRanges", "TranscriptDb"),
    function(query, subject, ...)
{
    .locateVariants(query, subject, ...)
})

.locateVariants <-
    function(query, subject, ..., cache=new.env(parent=emptyenv()))
{
    queryseq <- seqlevels(query)
    subseq <- seqlevels(subject)
    if (!any(queryseq %in% subseq))
        warning("none of seqlevels(query) match seqlevels(subject)")

    ## mask chromosomes not in query
    masks <- isActiveSeq(subject)
    on.exit(isActiveSeq(subject) <- masks)
    .setActiveSubjectSeq(query, subject)

    if (!exists(".__init__", cache, inherits=FALSE)) {
        cache[["tx"]] <- transcripts(subject,
            columns=c("tx_id", "gene_id"))
        cache[["cdsByTx"]] <- cdsBy(subject)
        cache[["fiveUTR"]] <- fiveUTRsByTranscript(subject)
        cache[["threeUTR"]] <- threeUTRsByTranscript(subject)
        cache[["txByGn"]] <- transcriptsBy(subject, "gene")
        cache[[".__init__"]] <- TRUE
    }
 
    map <- data.frame(
        txid=rep(names(cache[["cdsByTx"]]),
          elementLengths(cache[["cdsByTx"]])),
        cdsid=values(unlist(cache[["cdsByTx"]],
          use.names=FALSE))[["cds_id"]])

    ## ranges with width=0 :
    ## de-increment start to equal end value 
    if (any(insertion <- width(query) == 0))
        start(query)[insertion] <- start(query)[insertion] - 1

    cdsCO <- countOverlaps(query, cache[["cdsByTx"]], type="within")
    txFO <- findOverlaps(query, cache[["tx"]], type="within")
    txCO <- tabulate(queryHits(txFO), length(query))

    if (length(txFO) == 0) {
        mat1 <- DataFrame(queryID=integer(), location=.location(),
            txID=integer(), cdsID=character(), geneID=character(),
            precedesID=character(), followsID=character())
    } else {
        qhits <- queryHits(txFO)
        txID <- values(cache[["tx"]])["tx_id"][subjectHits(txFO),]
        cdsID <- map$cdsid[match(txID, map$txid)]
        geneID <- unlist(values(cache[["tx"]])[["gene_id"]][subjectHits(txFO)],
            use.names=FALSE)

        ## coding :
        coding <- cdsCO > 0 

        ## intron :
        intron <- txCO != 0 & cdsCO == 0

        ## UTRs :
        utr5 <-
            countOverlaps(query, cache[["fiveUTR"]], type="within") > 0
        utr3 <-
            countOverlaps(query, cache[["threeUTR"]], type="within") > 0

        location <- .location(length(qhits), "transcript_region")
        location[qhits %in% which(intron)] <- "intron"
        location[qhits %in% which(utr5)] <- "5'UTR"
        location[qhits %in% which(utr3)] <- "3'UTR"
        location[qhits %in% which(coding)] <- "coding"
        mat1 <- DataFrame(queryID=qhits, location, txID, cdsID, geneID,
            precedesID=NA_character_, followsID=NA_character_)
    }

    ## intergenic
    mat2 <- .intergenic(txCO, query, cache, map) 

    ans <- rbind(mat1, mat2)
    ans[order(ans$queryID), ]
}

.intergenic <- function(txCO, query, cache, map)
{
    intergenic <- txCO == 0
    if (all(intergenic == FALSE)) {
        DataFrame(queryID=integer(), location=.location(),
            txID=integer(), cdsID=integer(), geneID=character(), 
            precedesID=character(), followsID=character())
    } else {
        query <- query[intergenic]

        ## gene ranges
        rnglst <- range(cache[["txByGn"]])
        rng <- unlist(rnglst, use.names=FALSE)
        genes <- rep(names(rnglst), elementLengths(rnglst))
        ## query precedes subject; get index for following gene 
        pidx <- precede(query, rng)
        ## query follows subject; get index for preceding gene 
        fidx <- follow(query, rng)

        location <- .location(length(query))
        location[seqlevels(query) %in% seqlevels(rng)] <- "intergenic"

        DataFrame(queryID=which(intergenic), location=location, 
            txID=NA_integer_, cdsID=NA_integer_, geneID=NA_character_,
            precedesID=genes[pidx], followsID=genes[fidx])
    }
}
