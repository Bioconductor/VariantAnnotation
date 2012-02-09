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
            mat1 <- DataFrame(queryID=integer(), location=character(), 
                txID=integer(), geneID=CharacterList(), cdsID=integer())
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

            location <- rep("transcript_region", length(qhits))
            location[qhits %in% which(intron)] <- "intron"
            location[qhits %in% which(utr5)] <- "5'UTR"
            location[qhits %in% which(utr3)] <- "3'UTR"
            location[qhits %in% which(coding)] <- "coding"
            mat1 <- DataFrame(queryID=qhits, location, txID, geneID, cdsID)
        }

        ## intergenic :
        mat2 <- .intergenic(txCO, tx, query, subject, map) 

        ans <- rbind(mat1, mat2)
        ans$location <- factor(ans$location)
        ans[order(ans$queryID), ]
    }
)

.intergenic <- function(txCO, tx, query, subject, map, ...)
{
    if (!any(txCO == 0)) {
        DataFrame(queryID=integer(), location=character(), txID=integer(), 
            geneID=CharacterList(), cdsID=integer())
    } else {
        intergenic <- txCO == 0
        query <- query[txCO == 0]

        ## gene ranges
        txbygn <- transcriptsBy(subject, "gene")
        rnglst <- range(txbygn)
        rng <- unlist(rnglst, use.names=FALSE)
        genes <- rep(names(rnglst), elementLengths(rnglst))
        nidx <- nearest(query, rng)
        ## query precedes subject; get index for following gene 
        pidx <- precede(query, rng)
        ## query follows subject; get index for preceding gene 
        fidx <- follow(query, rng)

        geneid <- as.list(paste(genes[fidx], genes[pidx], sep=","))
        txid <- cdsid <- rep(NA_integer_, length(query))
        location <- rep(NA_character_, length(query))
        location[seqlevels(query) %in% seqlevels(rng)] <- "intergenic"

        DataFrame(queryID=which(intergenic), location=location, 
                  txID=txid, geneID=CharacterList(geneid), cdsID=cdsid) 
    }
}
