### =========================================================================
### locateVariants methods 
### =========================================================================

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

        ## ranges with width=0 :
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
            dat1 <- DataFrame(
                      queryHits=integer(), txID=character(), 
                      geneID=CharacterList(), Location=character())
        } else {
            queryHits <- queryHits(txFO)
            txID <- values(tx)["tx_id"][subjectHits(txFO),]
            geneID <- values(tx)[["gene_id"]][subjectHits(txFO)]
            if (any(elementLengths(geneID)) == 0)
                geneID[elementLengths(geneID) == 0] <- NA 

            ## coding :
            coding <- cdsCO > 0 

            ## intron :
            intron <- txCO != 0 & cdsCO == 0

            ## UTRs :
            fiveUTR <- fiveUTRsByTranscript(subject)
            threeUTR <- threeUTRsByTranscript(subject)
            utr5 <- countOverlaps(queryAdj, fiveUTR, type="within") > 0
            utr3 <- countOverlaps(queryAdj, threeUTR, type="within") > 0

            Location <- rep("transcript_region", length(queryHits))
            Location[queryHits %in% which(intron)] <- "intron"
            Location[queryHits %in% which(utr5)] <- "5'UTR"
            Location[queryHits %in% which(utr3)] <- "3'UTR"
            Location[queryHits %in% which(coding)] <- "coding"
            dat1 <- DataFrame(queryHits=queryHits, txID, geneID, Location)
        }

        ## intergenic :
        if (any(txCO == 0)) {
            intergenic <- txCO == 0
            ## tx that have gene IDs
            geneWidth <- elementLengths(values(tx)[["gene_id"]]) 
            if (any(geneWidth > 1))
                stop("assumption of one transcript per gene is not valid")
            txWithGeneID <- tx[geneWidth != 0]
            ## collapse to single range w/in gene ID
            grle <- Rle(unlist(values(txWithGeneID)[["gene_id"]], 
                use.names=FALSE))
            gfact <- rep(seq_len(nrun(grle)), runLength(grle)) 
            rngWithGeneID <- unlist(range(split(txWithGeneID, gfact)), use.names=FALSE)
            genes <- runValue(grle)

            ## locate range nearest to the variant 
            intvar <- queryAdj[txCO == 0]
            nst <- nearest(intvar, rngWithGeneID)
            if (any(is.na(nst))) {
                nonearestIdx <- which(intergenic == TRUE)[is.na(nst)] 
                nonearest <- rep(FALSE, length(intergenic))
                nonearest[nonearestIdx] <- TRUE
                intergenic[nonearestIdx] <- FALSE
                intvar <- queryAdj[intergenic == TRUE] 
                nearestIdx <- nearest(intvar, rngWithGeneID)
            } else {
                nearestIdx <- nearest(intvar, rngWithGeneID)
            }

            isPreceding <- (end(rngWithGeneID[nearestIdx]) - start(intvar)) < 0
            isFirst <- nearestIdx == 1
            isLast <- nearestIdx == length(rngWithGeneID) 

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
 
            flankGenes <- CharacterList(data.frame(rbind(unlist(p), unlist(f))))

            queryhits <- which(intergenic)
            txid <- rep(NA, length(which(intergenic)))
            geneid <- flankGenes
            location <- rep("intergenic", length(which(intergenic)))
            if (any(is.na(nst))) {
                queryhits <- c(queryhits, which(nonearest))
                txid <- c(txid, rep(NA, length(which(nonearest))))
                geneid <- c(flankGenes, 
                    CharacterList(as.list(rep(NA, length(which(nonearest)))))) 
                location <- c(location, rep(NA, length(which(nonearest))))
            }
            dat2 <- DataFrame(
                      queryHits=queryhits, txID=txid, geneID=geneid, 
                      Location=location)
        } else {
            dat2 <- DataFrame(
                      queryHits=integer(), txID=character(), 
                      geneID=CharacterList(), Location=character())
        }

        ans <- rbind(dat1, dat2)
        ans$Location <- factor(ans$Location)
        ans[order(ans$queryHits), ]
    }
)

