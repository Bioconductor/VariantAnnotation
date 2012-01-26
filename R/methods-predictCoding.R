### =========================================================================
### predictCoding methods 
### =========================================================================

setMethod("predictCoding",  signature(query="Ranges", subject="TranscriptDb", 
          seqSource="ANY", varAllele="DNAStringSet"),
    function(query, subject, seqSource, varAllele, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=subject, seqSource, 
            varAllele, ...) 
    }
)

setMethod("predictCoding", signature(query="VCF", subject="TranscriptDb", 
          seqSource="ANY", varAllele="missing"),
    function(query, subject, seqSource, varAllele, ...)
    {
        rd <- rowData(query) 
        alt <- values(alt(query))[["ALT"]]
        rdf <- rep(rd, elementLengths(alt))
        callGeneric(query=rdf, subject=subject, seqSource=seqSource, 
             varAllele=unlist(alt, use.names=FALSE), ...) 
    }
)

setMethod("predictCoding", signature(query="GRanges", subject="TranscriptDb", 
          seqSource="ANY", varAllele="DNAStringSet"),
    function(query, subject, seqSource, varAllele, ...)
    {
        stopifnot(length(varAllele) == length(query))

        ## check and adjust seqlevels
        txdb <- subject
        subject <- cdsBy(subject)
        qseq <- seqlevels(query)
        sseq <- seqlevels(subject)
        if (!any(qseq %in% sseq))
            warning("none of seqlevels(query) match seqlevels(subject)")
        subject <- keepSeqlevels(subject, query) 

        ## prevent findOverlaps circularity error
        circular <- .checkCircular(subject)
        if (!is.na(circular)) 
            stop("Please remove circular sequence ", circular, " from the ",
                 "query. Overlaps for type 'within' are not yet supported for ",
                 "circular sequences.")
 
        ## ranges with width=0 :
        ## de-increment start to equal end value 
        if (any(width(query) == 0)) {
            queryAdj <- query
            start(queryAdj[width(query) == 0]) <- 
                start(query)[width(query) == 0] - 1
        } else queryAdj <- query

        fo <- findOverlaps(queryAdj, subject, type = "within")
        if (length(fo) == 0)
            return(
                DataFrame(
                  queryID=character(0), consequence=character(0), 
                  refSeq=DNAStringSet(), varSeq=DNAStringSet(), 
                  refAA=AAStringSet(), varAA=AAStringSet(), 
                  txID=character(0), geneID=character(0), cdsID=character(0))) 

        ## map transcript indices to genome indices
        subject <- subject[unique(subjectHits(fo))]
        txLocal <- globalToLocal(queryAdj, subject)
        midx <- which(width(varAllele) == 0)
        if (length(midx) > 0) {
            warning("ranges with missing values in the varAllele column ",
                    "will be ignored")
            txLocal <- txLocal[!txLocal$globalInd %in% midx, ] 
        }
        xCoding <- query[txLocal$globalInd]
        xAllele <- varAllele[txLocal$globalInd]
 
        ## construct original sequences 
        originalWidth <- width(xCoding)
        codonStart <- (start(txLocal$local) - 1L) %/% 3L * 3L + 1L
        codonEnd <- codonStart + originalWidth %/% 3L * 3L + 2L
        txseqs <- getTranscriptSeqs(subject, seqSource)
        codons <- DNAStringSet(substring(txseqs[txLocal$rangesInd], 
            codonStart, codonEnd))

        ## construct variant sequences 
        varWidth <- width(xAllele)
        varPosition <- (start(txLocal$local) - 1L) %% 3L + 1L
        indels <- originalWidth == 0 | varWidth == 0 
        translateIdx <- abs(varWidth - originalWidth) %% 3 == 0 
        varSeq <- codons
        subseq(varSeq, start=varPosition, width=originalWidth) <- xAllele 

        ## results
        queryID <- txLocal$globalInd
        txID <- names(subject)[txLocal$rangesInd]
        fromSubject <-
            values(subject@unlistData)[txLocal$rangesInd,]
        if (any(names(fromSubject) %in% "cds_id"))
            cdsID <- fromSubject$cds_id
        else
            cdsID <- rep(NA, length(txID))
 
        txBygene <- transcriptsBy(txdb, "gene")
        map <- data.frame(
            geneid=rep(names(txBygene), elementLengths(txBygene)),
            txid=values(unlist(txBygene, use.names=FALSE))[["tx_id"]])
        geneID <- map$geneid[match(txID, map$txid)]

        refSeq <- codons
        varSeq <- varSeq
        refAA <- translate(codons)
        varAA <- AAStringSet(rep("", length(queryID))) 
        varAA[translateIdx] <- translate(varSeq[translateIdx])
 
        nonsynonymous <- as.character(refAA) != as.character(varAA) 
        consequence <- rep("synonymous", length(xCoding))
        consequence[nonsynonymous] <- "nonsynonymous" 
        consequence[!translateIdx] <- "frameshift" 
        consequence <- factor(consequence) 
 
        DataFrame(queryID, consequence, refSeq, varSeq, 
            refAA, varAA, txID, geneID, cdsID) 
    }
)

