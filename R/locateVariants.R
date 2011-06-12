locateVariants <- function(query, subject, ...)
{
t1 <- proc.time()
    chrom <- unique(seqlevels(query))
    #exons <- exons(subject, vals=list(exon_chrom=chrom),
    #   columns=c("exon_id", "tx_id", "gene_id"))
    cdsByTx <- cdsBy(subject)
    tx <- transcripts(subject, vals=list(exon_chrom=chrom),
       columns=c("exon_id", "tx_id", "gene_id"))
    #exFO <- findOverlaps(query, exons)
    cdsFO <- findOverlaps(query, cdsByTx, type="within")
    txFO <- findOverlaps(query, tx)
    #exCO <- tabulate(queryHits(exFO), length(query))
    cdsCO <- tabulate(queryHits(cdsFO), length(query))
    txCO <- tabulate(queryHits(txFO), length(query))
cat("set up", proc.time() - t1,"\n")

    ## exon :
    exonic <- cdsCO > 0 

t2 <- proc.time()
    ## intergenic :
    intergenic <- txCO == 0
    intvar <- query[txCO == 0]
    flankGenes <- lapply(seq_len(length(intvar)), 
        function(x, tx) {
            values(c(tx[precede(intvar[x,], tx)], 
            tx[follow(intvar[x,], tx)]))[["gene_id"]] 
        }, tx)
cat("intergenic",proc.time() - t2,"\n")

    ## intron :
    intronic <- txCO != 0 & cdsCO == 0

t4 <- proc.time()
    ## UTRs :
    fiveUTR <- fiveUTRsByTranscript(subject)
    threeUTR <- threeUTRsByTranscript(subject)
    utr5 <- countOverlaps(query, fiveUTR) > 0
    utr3 <- countOverlaps(query, threeUTR) > 0
cat("utr",proc.time() - t4,"\n")


t5 <- proc.time()
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
cat("geneID",proc.time() - t5,"\n")

    ## FIXME : report all matches or by priority?

    location <- character(length(query))
    location[exonic == TRUE] <- "exon"
    location[intronic == TRUE] <- "intron"
    location[intergenic == TRUE] <- "intergenic"
    location[utr5 == TRUE] <- "5'UTR"
    location[utr3 == TRUE] <- "3'UTR"

    values(query) <- append(values(query), DataFrame(location, txID, geneID))
    query
}


