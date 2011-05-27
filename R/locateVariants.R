locateVariants <- function(query, subject, orgPackage, ...)
{
    chrom <- unique(seqlevels(query))
    exons <- exons(subject, vals=list(exon_chrom=chrom),
       columns=c("exon_id", "tx_id", "gene_id"))
    tx <- transcripts(subject, vals=list(exon_chrom=chrom),
       columns=c("exon_id", "tx_id", "gene_id"))
    exFO <- findOverlaps(query, exons)
    txFO <- findOverlaps(query, tx)
    exCO <- tabulate(queryHits(exFO), length(query))
    txCO <- tabulate(queryHits(txFO), length(query))

    ## intergenic :
    intergenic <- txCO == 0
    intvar <- query[txCO == 0]
    flankGenes <- lapply(seq_len(length(intvar)), 
        function(x, tx) {
            values(c(tx[precede(intvar[x,], tx)], 
            tx[follow(intvar[x,], tx)]))[["gene_id"]] 
        }, tx)

    ## intronic :
    intronic <- txCO != 0 & exCO == 0

    ## UTRs :
    fiveUTR <- fiveUTRsByTranscript(subject)
    threeUTR <- threeUTRsByTranscript(subject)
    utr5 <- countOverlaps(query, fiveUTR) > 0
    utr3 <- countOverlaps(query, threeUTR) > 0

    ## exonic :
    exonic <- countOverlaps(query, cdsBy(subject), type="within") > 0 
    #exonic <- exCO > 0 & utr5 == FALSE & utr3 == FALSE

    ## gene metadata
    genes <- character(length(query)) 
    genes[txCO == 0] <- flankGenes
    genes[txCO != 0] <- split(values(tx)[["gene_id"]][subjectHits(txFO)], queryHits(txFO))
    genes <- lapply(genes, function(x) unique(as.numeric(x)))
    len <- unlist(lapply(genes, length), use.names=FALSE)

    ugenes <- unlist(genes, use.names=FALSE)
    x <- org.Hs.egSYMBOL
    mapped_symbols <- mappedkeys(x)
    xx <- as.list(x[mapped_symbols])
    gsymbol <- xx[names(xx) %in% ugenes]
    gsymbol <- gsymbol[match(ugenes, names(gsymbol))]
    #x <- org.Hs.egGENENAME
    #mapped_genes <- mappedkeys(x)
    #xx <- as.list(x[mapped_genes])
    #gname <- xx[names(xx) %in% ugenes]
    #gname <- gname[match(ugenes, names(gname))]
    #geneInfo <- unlist(paste(gsymbol, ": ", 
    #    gname[match(ugenes, names(gname))], sep=""), use.names=FALSE)
    genes <- CharacterList(split(gsymbol, rep(seq_len(length(len)), len)))

    ## FIXME : priority if >1 category=TRUE
    location <- character(length(query))
    location[exonic == TRUE] <- "exonic"
    location[intronic == TRUE] <- "intronic"
    location[intergenic == TRUE] <- "intergenic"
    location[utr5 == TRUE] <- "5'UTR"
    location[utr3 == TRUE] <- "3'UTR"

    values(query) <- append(values(query), DataFrame(location, genes))
    query
}


