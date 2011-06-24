locateVariants <- function(query, subject, ...)
{
    chrom <- unique(seqlevels(query))
    cdsByTx <- cdsBy(subject)
    tx <- transcripts(subject, vals=list(exon_chrom=chrom),
       columns=c("exon_id", "tx_id", "gene_id"))
    ## FIXME : type = 'within' or 'any' ?
    cdsFO <- findOverlaps(query, cdsByTx, type="within")
    txFO <- findOverlaps(query, tx)
    cdsCO <- tabulate(queryHits(cdsFO), length(query))
    txCO <- tabulate(queryHits(txFO), length(query))

    ## coding :
    coding <- cdsCO > 0 

    ## intergenic :
    ## FIXME : possible to not match a tx but be something other than intergenic? 
    if (any(txCO == 0)) {
        intergenic <- txCO == 0
        intvar <- query[txCO == 0]
        flankGenes <- lapply(seq_len(length(intvar)), 
            function(x, tx) {
                values(c(tx[precede(intvar[x,], tx)], 
                tx[follow(intvar[x,], tx)]))[["gene_id"]] 
            }, tx)
    } else flankGenes <- NA

    ## intron :
    intronic <- txCO != 0 & cdsCO == 0

    ## UTRs :
    fiveUTR <- fiveUTRsByTranscript(subject)
    threeUTR <- threeUTRsByTranscript(subject)
    utr5 <- countOverlaps(query, fiveUTR) > 0
    utr3 <- countOverlaps(query, threeUTR) > 0

    ## geneID 
    genes <- character(length(query)) 
    transcripts <- character(length(query)) 
    genes[txCO == 0] <- flankGenes
    genes[txCO != 0] <- split(values(tx)[["gene_id"]][subjectHits(txFO)], queryHits(txFO))
    transcripts[txCO != 0] <- split(values(tx)[["tx_id"]][subjectHits(txFO)], queryHits(txFO))
    genes <- lapply(genes, function(x) unique(as.numeric(x)))
    transcripts <- lapply(transcripts, function(x) unique(as.numeric(x)))
    geneID <- CharacterList(genes)
    txID <- CharacterList(transcripts)

    ## FIXME : report all matches or only priority?
    Location <- rep("unknown", length(query))
    Location[intergenic == TRUE] <- "intergenic"
    Location[intronic == TRUE] <- "intron"
    Location[utr5 == TRUE] <- "5'UTR"
    Location[utr3 == TRUE] <- "3'UTR"
    Location[coding == TRUE] <- "coding"

    values(query) <- append(values(query), DataFrame(Location, txID, geneID))
    query
}


