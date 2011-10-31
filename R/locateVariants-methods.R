### =========================================================================
### locateVariants methods 
### =========================================================================

setMethod("locateVariants", c("GRanges", "TranscriptDb"),
    function(query, subject, ...)
    {
        chrom <- seqlevels(query)
        isActiveSeq(subject)[!names(isActiveSeq(subject)) %in% chrom] <- FALSE 
        tx <- transcripts(subject, columns=c("exon_id", "tx_id", "gene_id"))
        cdsByTx <- cdsBy(subject)

        ## findOverlaps won't find negative widths
        ## adjust query with width=0 :
        ## de-increment start to equal end value 
        if (any(width(query) == 0)) {
            insertion <- width(query) == 0
            queryAdj <- query
            start(queryAdj[insertion]) <- 
                start(query)[insertion] - 1
        } else queryAdj <- query

        cdsFO <- findOverlaps(queryAdj, cdsByTx, type="within")
        txFO <- findOverlaps(queryAdj, tx, type="within")
        cdsCO <- tabulate(queryHits(cdsFO), length(query))
        txCO <- tabulate(queryHits(txFO), length(query))
        if (sum(txCO) == 0) {
            return(
              DataFrame(queryHits=seq_len(length(query)),
                        txID=character(length(query)),
                        geneID=CharacterList(character(length(query))),
                        Location=factor(rep("unknown", length(query))))
            )
        }

        queryHits <- queryHits(txFO)
        txID <- names(cdsByTx)[subjectHits(txFO)]
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

        ## intergenic :
        if (any(txCO == 0)) {
            intergenic <- txCO == 0
            intvar <- queryAdj[txCO == 0]
            geneWidth <- width(values(tx)[["gene_id"]]@partitioning) 
            txWithGeneID <- tx[geneWidth != 0]
            genes <- values(txWithGeneID)[["gene_id"]]

            nearestIdx <- nearest(intvar, txWithGeneID)
            isPreceding <- (end(txWithGeneID[nearestIdx]) - start(intvar)) < 0
            isFirst <- nearestIdx == 1
            isLast <- nearestIdx == length(txWithGeneID) 

            ## nearest is preceding
            precedeIdx <- nearestIdx
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

            ## assuming a transcript can fall in only 1 gene 
            flankGenes <- CharacterList(data.frame(rbind(unlist(p), unlist(f))))
            dat2 <- 
              DataFrame(queryHits=which(intergenic), 
                        txID=rep(NA, length(which(intergenic))), 
                        geneID=flankGenes, 
                        Location=rep("intergenic", length(which(intergenic))))
        } else {
            dat2 <- 
              DataFrame(queryHits=integer(), txID=character(), 
                        geneID=CharacterList(), Location=character())
        }

        Location <- rep("unknown", length(queryHits))
        Location[queryHits %in% which(intron)] <- "intron"
        Location[queryHits %in% which(utr5)] <- "5'UTR"
        Location[queryHits %in% which(utr3)] <- "3'UTR"
        Location[queryHits %in% which(coding)] <- "coding"
        dat1 <- DataFrame(queryHits=queryHits, txID, geneID, Location)
        ans <- rbind(dat1, dat2)
        ans$Location <- factor(ans$Location)
        ans[order(ans$queryHits), ]
    }
)


