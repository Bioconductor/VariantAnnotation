setMethod("predictCoding",  c("Ranges", "TranscriptDb"),
    function(query, subject, variantAllele, seqSource, ...)
    {
        ## FIXME : pass on only subj that hits query
        cdsByTx <- cdsBy(subject)
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=cdsByTx, variantAllele=variantAllele,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding",  c("GRanges", "TranscriptDb"),
    function(query, subject, variantAllele, seqSource, ...)
    {
        ## FIXME : pass on only subj that hits query
        cdsByTx <- cdsBy(subject)
        callGeneric(query=query, subject=cdsByTx, variantAllele=variantAllele,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding", c("Ranges", "GRangesList"),
    function(query, subject, variantAllele, seqSource, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=subject, variantAllele=variantAllele,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding", c("GRanges", "GRangesList"),
    function(query, subject, variantAllele, seqSource, ...)
    {
        ## FIXME : should query include only variants in exons?
        ##         findOverlaps is repeated 3x
        txSeqs <- getTranscriptSeqs(subject, seqSource)
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
            as.character(values(xCoding)[[variantAllele]])

        values(xCoding)[colnames(values(subject))] <-
            values(subject)[txLocal$ranges.ind,]
        values(xCoding)$txID <- names(subject)[txLocal$ranges.ind]
        values(xCoding)$refAA <- translate(DNAStringSet(codons))
        values(xCoding)$obsAA <- translate(DNAStringSet(varCodons))
        xCoding
    }
)

