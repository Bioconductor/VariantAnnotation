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

        queryseq <- seqlevels(query)
        subseq <- seqlevels(subject)
        if (!any(queryseq %in% subseq))
            warning("none of seqlevels(query) match seqlevels(subject)")

        ## mask chromosomes not in query
        masks <- isActiveSeq(subject)
        isActiveSeq(subject)[!names(isActiveSeq(subject)) %in% queryseq] <- FALSE
        on.exit(isActiveSeq(subject) <- masks)

        ## prevent findOverlaps circularity error
        circular <- .checkCircular(subject)
        if (any(circular %in% isActiveSeq(subject))) 
            stop("Please remove circular sequence ", circular, " from the ",
                 "query. Overlaps for type 'within' are not yet supported for ",
                 "circular sequences.")

        cdsByTx <- cdsBy(subject)
        ## ranges with width=0 :
        ## de-increment start to equal end value 
        if (any(width(query) == 0)) {
            queryAdj <- query
            start(queryAdj[width(query) == 0]) <- 
                start(query)[width(query) == 0] - 1
        } else queryAdj <- query

        fo <- findOverlaps(queryAdj, cdsByTx, type = "within")
        if (length(fo) == 0)
            return(
                DataFrame(
                  queryID=integer(), consequence=character(), 
                  refSeq=DNAStringSet(), varSeq=DNAStringSet(), 
                  refAA=AAStringSet(), varAA=AAStringSet(), 
                  txID=integer(), geneID=character(), cdsID=integer())) 

        ## map transcript indices to genome indices
        cdsByTx <- cdsByTx[unique(subjectHits(fo))]
        txLocal <- globalToLocal(queryAdj, cdsByTx)
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
        txseqs <- getTranscriptSeqs(cdsByTx, seqSource)
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
        txID <- names(cdsByTx)[txLocal$rangesInd]
        fromSubject <-
            values(cdsByTx@unlistData)[txLocal$rangesInd,]
        if (any(names(fromSubject) %in% "cds_id"))
            cdsID <- fromSubject$cds_id
        else
            cdsID <- rep(NA, length(txID))
 
        txByGene <- transcriptsBy(subject, "gene")
        map <- data.frame(
            geneid=rep(names(txByGene), elementLengths(txByGene)),
            txid=values(unlist(txByGene, use.names=FALSE))[["tx_id"]])
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

