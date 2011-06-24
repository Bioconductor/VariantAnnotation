setMethod("predictCoding",  c("Ranges", "TranscriptDb"),
    function(query, subject, alleles, seqSource, ...)
    {
        cdsByTx <- cdsBy(subject)
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=cdsByTx, alleles=alleles,
            seqSource=seqSource, ...) 
    }
)

setMethod("predictCoding",  c("GRanges", "TranscriptDb"),
    function(query, subject, alleles, seqSource, ...)
    {
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
        ## FIXME : findOverlaps is done here, globalToLocal and locateVariants
        fo <- findOverlaps(query, subject, type = "within")
        subject <- subject[unique(subjectHits(fo))]
        txSeqs <- getTranscriptSeqs(subject, seqSource)
        txLocal <- globalToLocal(query, subject)
        xCoding <- query[txLocal$global.ind]
        splitAlleles <- do.call(rbind, strsplit(values(xCoding)[[alleles]], "/"))
        refAlleles <- splitAlleles[,1] 
        varAlleles <- splitAlleles[,2] 

        ## FIXME : check original sequence from BSgenome/fasta
        ##         against user provided refAllele?
 
        ## original codons
        originalWidth <- width(xCoding)
        codonStart <- (start(txLocal$local) - 1L) %/% 3L * 3L + 1L
        codonEnd <- codonStart + 
            originalWidth %/% 3L * 3L + 2L
        codons <- tolower(substring(txSeqs[txLocal$ranges.ind], codonStart,
            codonEnd))
	#substring(codons, varPosition, varPosition) <- 
	#    toupper(substring(codons, varPosition, varPosition)) 

        ## variant codons
        varWidth <- nchar(varAlleles)
        varPosition <- (start(txLocal$local) - 1L) %% 3L + 1L
        insertions <- refAlleles == "-" 
        deletions <- varAlleles == "-" 
        indels <- insertions | deletions
        allowedSubstitutions <- !indels &
            abs(varWidth - originalWidth) %% 3L == 0
        allowedIndels <- indels  &
            (abs(varWidth - originalWidth) + 1L) %% 3L == 0
        translateIdx <- allowedSubstitutions | allowedIndels 
        varCodons <- codons


	## FIXME : simplify
        ## substitutions, including snvs
        first <- substr(varCodons[!indels], 1,
            varPosition[!indels] - 1L)
        last <- substr(varCodons[!indels], varPosition[!indels] + 1L,
            nchar(varCodons[!indels]))
        varCodons[!indels] <- paste(first, varAlleles[!indels], last, sep="")
        ## insertions
        first <- substr(varCodons[insertions], varPosition[insertions] - 1L,
            varPosition[insertions] - 1L)
        last <- substr(varCodons[insertions], varPosition[insertions],
            nchar(varCodons[insertions]))
        varCodons[insertions] <- paste(first, varAlleles[insertions], last, sep="")
        ## deletions
        first <- substr(varCodons[deletions], 1, varPosition[deletions] - 1L)
        last <- substr(varCodons[deletions], varPosition[deletions] +
            originalWidth[deletions], nchar(varCodons[deletions]))
        varCodons[deletions] <- paste(first, last, sep="")
	

        ## results
        values(xCoding)[colnames(values(subject@unlistData))] <-
            values(subject@unlistData)[txLocal$ranges.ind,]
        values(xCoding)$tx_ID <- names(subject)[txLocal$ranges.ind]

        CodonChange <- paste(codons, "/",
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

