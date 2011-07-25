## TODO : (1) Currently overlaps all type=within.
##            type=any and list all matches?
##            Note that predictCoding uses type=any
##        (2) Simplify flanking genes search for intergenic
##        (3) need findOverlaps() method for decimal location 
##            (ie, 3.4 is insertion location btn 3 and 4)

setMethod("locateVariants", c("GRanges", "TranscriptDb"),
    function(query, subject, ...)
    {
        chrom <- seqlevels(query)
        #isActiveSeq(subject)[!names(isActiveSeq(subject)) %in% chrom] <- FALSE 
        #tx <- transcripts(subject, columns=c("exon_id", "tx_id", "gene_id"))
        cdsByTx <- cdsBy(subject)
        tx <- transcripts(subject, vals=list(exon_chrom=chrom),
           columns=c("exon_id", "tx_id", "gene_id"))


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

        Location <- rep("unknown", length(queryHits))
        Location[queryHits %in% which(intron)] <- "intron"
        Location[queryHits %in% which(utr5)] <- "5'UTR"
        Location[queryHits %in% which(utr3)] <- "3'UTR"
        Location[queryHits %in% which(coding)] <- "coding"

        dat1 <- DataFrame(queryHits=queryHits, txID, geneID, Location)

        ## intergenic :
        ## FIXME : need to simplify
        if (any(txCO == 0)) {
            intergenic <- txCO == 0
            intvar <- query[txCO == 0]
            flankGenes <- lapply(seq_len(length(intvar)), 
                function(x, tx, intvar) {
                    pre <- fol <- intvar[x,]
                    pregene <- values(tx[precede(pre, tx)])["gene_id"] 
                    while (length(pregene$gene_id[[1]]) == 0 && 
                           (typeof(pregene$gene_id[[1]]) == "character")){
                            pre <- precede(pre, tx)
                            pregene <- values(tx[precede(pre, tx)])["gene_id"] 
                    }

                    folgene <- values(tx[follow(fol, tx)])["gene_id"] 
                    while (length(folgene$gene_id[[1]]) == 0 && 
                           (typeof(folgene$gene_id[[1]]) == "character")){
                            fol <- tx[follow(fol, tx)]
                            folgene <- values(tx[follow(fol, tx)])["gene_id"] 
                    }
                c(pregene$gene_id[[1]], folgene$gene_id[[1]])
                }, tx, intvar)
        } else flankGenes <- NA

       dat2 <- DataFrame(queryHits=which(intergenic), 
           txID=rep(NA, length(which(intergenic))), 
           geneID=CharacterList(flankGenes), 
           Location=rep("intergenic", length(which(intergenic))))

       ans <- rbind(dat1, dat2)
       ans[order(ans$queryHits), ]
    }
)


