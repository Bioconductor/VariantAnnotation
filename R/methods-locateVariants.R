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
        isActiveSeq(subject)[!names(isActiveSeq(subject)) %in% queryseq] <- FALSE
        on.exit(isActiveSeq(subject) <- masks)

        ## prevent findOverlaps circularity error
        circular <- .checkCircular(subject)
        if (any(circular %in% isActiveSeq(subject))) 
            stop("Please remove circular sequence ", circular, " from the ",
                 "query. Overlaps for type 'within' are not yet supported for ",
                 "circular sequences.")

        tx <- transcripts(subject, columns=c("exon_id", "tx_id", "gene_id"))
        cdsByTx <- cdsBy(subject)
        map <- data.frame(
            txid=rep(names(cdsByTx), elementLengths(cdsByTx)),
            cdsid=values(unlist(cdsByTx, use.names=FALSE))[["cds_id"]])

        ## ranges with width=0 :
        ## de-increment start to equal end value 
        if (any(width(query) == 0)) {
            insertion <- width(query) == 0
            queryAdj <- query
            start(queryAdj[insertion]) <- start(query)[insertion] - 1
        } else queryAdj <- query

        cdsFO <- findOverlaps(queryAdj, cdsByTx, type="within")
        cdsCO <- tabulate(queryHits(cdsFO), length(queryAdj))
        txFO <- findOverlaps(queryAdj, tx, type="within")
        txCO <- tabulate(queryHits(txFO), length(queryAdj))

        if (sum(txCO) == 0) {
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
            threeUTR <- threeUTRsByTranscript(subject)
            utr5 <- countOverlaps(queryAdj, fiveUTR, type="within") > 0
            utr3 <- countOverlaps(queryAdj, threeUTR, type="within") > 0

            location <- rep("transcript_region", length(qhits))
            location[qhits %in% which(intron)] <- "intron"
            location[qhits %in% which(utr5)] <- "5'UTR"
            location[qhits %in% which(utr3)] <- "3'UTR"
            location[qhits %in% which(coding)] <- "coding"
            mat1 <- DataFrame(queryID=qhits, location, txID, geneID, cdsID)
        }

        ## intergenic :
        mat2 <- .intergenic(txCO, tx, queryAdj, map) 

        ans <- rbind(mat1, mat2)
        ans$location <- factor(ans$location)
        ans[order(ans$queryID), ]
    }
)

.intergenic <- function(txCO, tx, queryAdj, map, ...)
{
    if (!any(txCO == 0)) {
        DataFrame(queryID=integer(), location=character(), txID=integer(), 
            geneID=CharacterList(), cdsID=integer())
    } else {
        intergenic <- txCO == 0
        ## tx that have gene IDs
        geneWidth <- elementLengths(values(tx)[["gene_id"]]) 
        if (any(geneWidth > 1))
            stop("assumption of one gene per transcript is not valid")
        txWithGeneID <- tx[geneWidth != 0]
        ## collapse to single range w/in gene ID
        grle <- Rle(unlist(values(txWithGeneID)[["gene_id"]], 
            use.names=FALSE))
        gfact <- rep(seq_len(nrun(grle)), runLength(grle)) 
        rngWithGeneID <- unlist(range(split(txWithGeneID, gfact)), 
            use.names=FALSE)
        genes <- runValue(grle)

        ## locate nearest range 
        intvar <- queryAdj[txCO == 0]
        nidx <- nearest(intvar, rngWithGeneID)
        nnidx <- logical(length(txCO))
        if (any(is.na(nidx))) {
            nnidx[which(txCO == 0)] <- is.na(nidx)
            intergenic[nnidx] <- FALSE
            intvar <- intvar[!is.na(nidx)]
            nidx <- na.omit(nidx)
        } 
        isPreceding <- (end(rngWithGeneID[nidx]) - start(intvar)) < 0
        isFirst <- nidx == 1
        isLast <- nidx == length(rngWithGeneID) 
        ## nearest is preceding
        precedeIdx <- nidx
        followIdx <- precedeIdx + 1
        ## nearest is following 
        precedeIdx[!isPreceding & !isFirst] <- 
            precedeIdx[!isPreceding & !isFirst] - 1 
        followIdx[!isPreceding & !isFirst] <- 
            followIdx[!isPreceding & !isFirst] - 1 

        ## edge cases where nearest is first or last 
        ## (ie, one flanking gene)
        p <- genes[precedeIdx]
        f <- genes[followIdx]
        if (any(isPreceding & isLast)) 
            f[isPreceding & isLast] <- "no following" 
        if (any(!isPreceding & isFirst)) 
            p[!isPreceding & isFirst] <- "no preceding"
        flankGenes <- CharacterList(data.frame(rbind(unlist(p), unlist(f))))

        qhits <- which(intergenic)
        txid <- cdsid <- rep(NA_integer_, length(qhits))
        geneid <- flankGenes
        location <- rep("intergenic", length(qhits))
        res <- DataFrame(queryID=qhits, location=location, txID=txid,
            geneID=geneid, cdsID=cdsid) 
        if (any(nnidx)) { 
            noNearestDF <- DataFrame(queryID=which(nnidx), location=NA_character_, 
                txID=NA_integer_, geneID=CharacterList(NA_character_), cdsID=NA_integer_)
            DF <- rbind(res, noNearestDF)
            DF[order(DF$queryID), ] 
        } else {
            res
        }
    }
}
