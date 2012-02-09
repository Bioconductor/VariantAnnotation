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
        queryseq <- seqlevels(query)
        subseq <- seqlevels(subject)
        if (!any(queryseq %in% subseq))
            warning("none of seqlevels(query) match seqlevels(subject)")

        ## mask chromosomes not in query
        masks <- isActiveSeq(subject)
        on.exit(isActiveSeq(subject) <- masks)
        .setActiveSubjectSeq(query, subject)

        tx <- transcripts(subject, columns=c("exon_id", "tx_id", "gene_id"))
        cdsByTx <- cdsBy(subject)
        map <- data.frame(
            txid=rep(names(cdsByTx), elementLengths(cdsByTx)),
            cdsid=values(unlist(cdsByTx, use.names=FALSE))[["cds_id"]])

        ## ranges with width=0 :
        ## de-increment start to equal end value 
        if (any(insertion <- width(query) == 0))
            start(query)[insertion] <- start(query)[insertion] - 1

        cdsCO <- countOverlaps(query, cdsByTx, type="within")
        txFO <- findOverlaps(query, tx, type="within")
        txCO <- tabulate(queryHits(txFO), length(query))

        if (length(txFO) == 0) {
            mat1 <- DataFrame(queryID=integer(), location=.location(),
                txID=integer(), geneID=character(),
                precedesID=character(), followsID=character(),
                cdsID=integer())
        } else {
            qhits <- queryHits(txFO)
            txID <- values(tx)["tx_id"][subjectHits(txFO),]
            cdsID <- map$cdsid[match(txID, map$txid)]
            geneID <- values(tx)[["gene_id"]][subjectHits(txFO)]

            ## coding :
            coding <- cdsCO > 0 

            ## intron :
            intron <- txCO != 0 & cdsCO == 0

            ## UTRs :
            fiveUTR <- fiveUTRsByTranscript(subject)
            utr5 <- countOverlaps(query, fiveUTR, type="within") > 0
            threeUTR <- threeUTRsByTranscript(subject)
            utr3 <- countOverlaps(query, threeUTR, type="within") > 0

            location <- .location(length(qhits), "transcript_region")
            location[qhits %in% which(intron)] <- "intron"
            location[qhits %in% which(utr5)] <- "5'UTR"
            location[qhits %in% which(utr3)] <- "3'UTR"
            location[qhits %in% which(coding)] <- "coding"
            mat1 <- DataFrame(queryID=qhits, location, txID, geneID,
                precedesID=NA_character_, followsID=NA_character_, cdsID)
        }

        ## intergenic :
        mat2 <- .intergenic(txCO, tx, query, subject, map) 

        ans <- rbind(mat1, mat2)
        ans[order(ans$queryID), ]
    }
)

.intergenic <- function(txCO, tx, query, subject, map, ...)
{
    intergenic <- txCO == 0
    if (all(intergenic)) {
        DataFrame(queryID=integer(), location=.location(),
            txID=integer(), geneID=character(), precedesID=character(),
            followsID=character(), cdsID=integer())
    } else {
        query <- query[intergenic]

        ## gene ranges
        txbygn <- transcriptsBy(subject, "gene")
        rnglst <- range(txbygn)
        rng <- unlist(rnglst, use.names=FALSE)
        genes <- rep(names(rnglst), elementLengths(rnglst))
        ## query precedes subject; get index for following gene 
        pidx <- precede(query, rng)
        ## query follows subject; get index for preceding gene 
        fidx <- follow(query, rng)

        location <- .location(length(query))
        location[seqlevels(query) %in% seqlevels(rng)] <- "intergenic"

        DataFrame(queryID=which(intergenic), location=location, 
            txID=NA_integer_, geneID=NA_character_,
            precedesID=genes[pidx], followsID=genes[fidx],
            cdsID=NA_integer_) 
    }
}
