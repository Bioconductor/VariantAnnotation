setMethod("predictCoding",  c("Ranges", "TranscriptDb"),
    function(query, subject, alleles, seqSource, ...)
    {
        ## FIXME : pass on only subj that hits query
        cdsByTx <- cdsBy(subject)
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=cdsByTx, alleles=alleles,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding",  c("GRanges", "TranscriptDb"),
    function(query, subject, alleles, seqSource, ...)
    {
        ## FIXME : pass on only subj that hits query
        cdsByTx <- cdsBy(subject)
        callGeneric(query=query, subject=cdsByTx, alleles=alleles,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding", c("Ranges", "GRangesList"),
    function(query, subject, alleles, seqSource, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=subject, alleles=alleles,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding", c("GRanges", "GRangesList"),
    function(query, subject, alleles, seqSource, ...)
    {
        txSeqs <- getTranscriptSeqs(subject, seqSource)
        ## FIXME : query include only variants in exons?
        ##         findOverlaps is repeated 3x
        txLocal <- globalToLocal(query, subject)
        xCoding <- query[txLocal$global.ind]
        splitAlleles <- do.call(rbind, strsplit(values(xCoding)[[alleles]], "/"))
        refAlleles <- splitAlleles[,1] 
        varAlleles <- splitAlleles[,2] 

        ## FIXME : check original sequence from BSgenome/fasta
        ##         against user provided refAllele?
        ##         check observed allele is same length as range
 
        ## original codons
        originalWidth <- width(xCoding)
        codonStart <- (start(txLocal$local) - 1L) %/% 3L * 3L + 1L
        codonEnd <- codonStart + 
            originalWidth %/% 3L * 3L + 2L
        codons <- substring(txSeqs[txLocal$ranges.ind], codonStart, codonEnd)

        ## variant codons
        ## FIXME : treat indels like substitutions?
        varWidth <- nchar(varAlleles)
        varPosition <- (start(txLocal$local) - 1L) %% 3L + 1L
        substitutionIdx <- (varWidth - originalWidth) %% 3L == 0
        indelIdx <- refAlleles == "-" | varAlleles == "-" 
        translateIdx <- substitutionIdx & !indelIdx 
        varCodons <- tolower(codons)
	substring(varCodons, varPosition, (varPosition + varWidth - 1L)) <- 
            toupper(as.character(varAlleles))

        ## results
        values(xCoding)[colnames(values(subject@unlistData))] <-
            values(subject@unlistData)[txLocal$ranges.ind,]
        values(xCoding)$tx_ID <- names(subject)[txLocal$ranges.ind]

        CodonChange <- paste(tolower(codons), "/",
            varCodons, sep="") 
        CodonChange[!translateIdx] <- NA
        values(xCoding)$CodonChange <- CodonChange 

        refAA <- translate(DNAStringSet(codons))
        varAA <- NA 
        varAA[translateIdx] <-
            as.character(translate(DNAStringSet(varCodons)[translateIdx]))
        AAChange <- paste(refAA, "/", varAA, sep="")
        AAChange[!translateIdx] <- NA 
        values(xCoding)$AAChange <- AAChange

        ## FIXME : better check for equality
        nonsynonymous <- as.character(refAA) != varAA 
        Consequence <- rep("synonymous_coding", length(xCoding))
        Consequence[nonsynonymous] <- "nonsynonymous_coding" 
        Consequence[!translateIdx] <- "frameshift_coding" 
        values(xCoding)$Consequence <- Consequence 
        
        xCoding
    }
)

