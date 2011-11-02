### =========================================================================
### predictCoding methods 
### =========================================================================

setMethod("predictCoding",  signature(query="Ranges", subject="TranscriptDb", 
          seqSource="ANY", varAllele="character"),
    function(query, subject, seqSource, varAllele, ...)
    {
        cdsByTx <- cdsBy(subject)
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=cdsByTx, seqSource, 
            varAllele, ...) 
    }
)

setMethod("predictCoding",  signature(query="GRanges", subject="TranscriptDb", 
          seqSource="ANY", varAllele="character"),
    function(query, subject, seqSource, varAllele, ...)
    {
        cdsByTx <- cdsBy(subject)
        callGeneric(query=query, subject=cdsByTx, seqSource, 
            varAllele, ...) 
    }
)

setMethod("predictCoding", signature(query="Ranges", subject="GRangesList", 
          seqSource="ANY", varAllele="character"),
    function(query, subject, seqSource, varAllele, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=subject, seqSource=seqSource, 
             varAllele=varAllele, ...) 
    }
)

setMethod("predictCoding", signature(query="GRanges", subject="GRangesList", 
          seqSource="ANY", varAllele="character"),
    function(query, subject, seqSource, varAllele, ...)
    {
        ## findOverlaps won't find negative widths
        ## adjust query with width=0 :
        ## de-increment start to equal end value 
        if (any(width(query) == 0)) {
            queryAdj <- query
            start(queryAdj[width(query) == 0]) <- 
                start(query)[width(query) == 0] - 1
        } else queryAdj <- query
 
        if (is.null(values(query)[[varAllele]]))
            stop("varAllele column not present in query")
 
        fo <- findOverlaps(queryAdj, subject, type = "within")
        if (length(fo) == 0)
        return(DataFrame(queryHits=character(0), txID=character(0),
               refSeq=DNAStringSet(), varSeq=DNAStringSet(), 
               refAA=AAStringSet(), varAA=AAStringSet(), 
               Consequence=character(0))) 

        subject <- subject[unique(subjectHits(fo))]
        txSeqs <- getTranscriptSeqs(subject, seqSource)
        txLocal <- globalToLocal(queryAdj, subject)
        xCoding <- query[txLocal$globalInd]
 
        ## original sequences 
        originalWidth <- width(xCoding)
        codonStart <- (start(txLocal$local) - 1L) %/% 3L * 3L + 1L
        codonEnd <- codonStart + 
            originalWidth %/% 3L * 3L + 2L
        codons <- DNAStringSet(substring(txSeqs[txLocal$rangesInd], 
            codonStart, codonEnd))

        ## variant sequences 
        varWidth <- width(values(xCoding)[[varAllele]])
        varPosition <- (start(txLocal$local) - 1L) %% 3L + 1L
        indels <- originalWidth == 0 | varWidth == 0 
        translateIdx <- abs(varWidth - originalWidth) %% 3 == 0 
        varSeq <- codons
        subseq(varSeq, start=varPosition, width=originalWidth) <- 
            values(xCoding)[[varAllele]]

        ## results
        queryHits <- txLocal$globalInd
        txID <- names(subject)[txLocal$rangesInd]
        fromSubject <-
            values(subject@unlistData)[txLocal$rangesInd,]
        refSeq <- codons
        varSeq <- varSeq
        refAA <- translate(codons)
        varAA <- AAStringSet(rep("", length(queryHits))) 
        varAA[translateIdx] <- translate(varSeq[translateIdx])
 
        nonsynonymous <- as.character(refAA) != as.character(varAA) 
        Consequence <- rep("synonymous", length(xCoding))
        Consequence[nonsynonymous] <- "nonsynonymous" 
        Consequence[!translateIdx] <- "frameshift" 
        Consequence <- factor(Consequence) 
 
        DataFrame(queryHits, txID, refSeq, varSeq, 
            refAA, varAA, Consequence, fromSubject) 
    }
)

