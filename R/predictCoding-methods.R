setMethod("predictCoding",  c("Ranges", "TranscriptDb"),
    function(query, subject, variantColumn, seqSource, ...)
    {
        ## FIXME : pass on only subj that hits query
        cdsByTx <- cdsBy(subject)
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=cdsByTx, variantColumn=variantColumn,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding",  c("GRanges", "TranscriptDb"),
    function(query, subject, variantColumn, seqSource, ...)
    {
        ## FIXME : pass on only subj that hits query
        cdsByTx <- cdsBy(subject)
        callGeneric(query=query, subject=cdsByTx, variantColumn=variantColumn,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding", c("Ranges", "GRangesList"),
    function(query, subject, variantColumn, seqSource, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=subject, variantColumn=variantColumn,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding", c("GRanges", "GRangesList"),
    function(query, subject, variantColumn, seqSource, ...)
    {
        ## FIXME : should query include only variants in exons?
        ##         findOverlaps is repeated 3x
        txSeqs <- getTranscriptSeqs(seqSource, subject)
        txLocal <- globalToLocal(query, subject)
        xCoding <- query[txLocal$global.ind]
        ## extract codon sequences
        codonStart <- (start(txLocal$local) - 1L) %/% 3L * 3L + 1L
        codons <- substring(txSeqs[txLocal$ranges.ind], codonStart, 
            codonStart + 2L)
        codonPos <- (start(txLocal$local) - 1L) %% 3L + 1L

        ## assuming variant with single base pair change (ie, snv)
        varCodons <- codons
        substring(varCodons, codonPos, codonPos) <- 
            as.character(values(xCoding)[[variantColumn]])

        values(xCoding)[colnames(values(subject))] <-
            values(subject)[txLocal$ranges.ind,]
        values(xCoding)$txID <- names(subject)[txLocal$ranges.ind]
        values(xCoding)$refAA <- translate(DNAStringSet(codons))
        values(xCoding)$obsAA <- translate(DNAStringSet(varCodons))
        xCoding
    }
)

