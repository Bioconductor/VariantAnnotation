setMethod("predictCoding",  c("Ranges", "TranscriptDb"),
    function(query, subject, variantColumn, BSgenomeOrganism, ...)
    {
        cdsByTx <- cdsBy(subject)
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=cdsByTx, variantColumn=variantColumn,
            BSgenomeOrganism=BSgenomeOrganism, ...) 
    }
)

setMethod("predictCoding",  c("Ranges", character(0)),
    function(query, subject, variantColumn, BSgenomeOrganism, ...)
    {
        ## FIXME : bam method
        #cdsByTx <- cdsBy(subject)
        #callGeneric(query=query, subject=cdsByTx, variantColumn=variantColumn,
        #    BSgenomeOrganism=BSgenomeOrganism, ...) 
    }
)

setMethod("predictCoding",  c("GRanges", "TranscriptDb"),
    function(query, subject, variantColumn, BSgenomeOrganism, ...)
    {
        cdsByTx <- cdsBy(subject)
        callGeneric(query=query, subject=cdsByTx, variantColumn=variantColumn,
            BSgenomeOrganism=BSgenomeOrganism, ...) 
    }
)

setMethod("predictCoding",  c("GRanges", character(0)),
    function(query, subject, variantColumn, BSgenomeOrganism, ...)
    {
        ## FIXME : bam method
        #cdsByTx <- cdsBy(subject)
        #callGeneric(query=query, subject=cdsByTx, variantColumn=variantColumn,
        #    BSgenomeOrganism=BSgenomeOrganism, ...) 
    }
)

setMethod("predictCoding", c("Ranges", "GRangesList"),
    function(query, subject, variantColumn, BSgenomeOrganism, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=subject, variantColumn=variantColumn,
            BSgenomeOrganism=BSgenomeOrganism, ...) 
    }
)

setMethod("predictCoding", c("GRanges", "GRangesList"),
    function(query, subject, variantColumn, BSgenomeOrganism, ...)
    {
        txSeqs <- getTranscriptSeqs(subject, BSgenomeOrganism) 
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

