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
        res <- callGeneric(query=rdf, subject=subject, seqSource=seqSource, 
             varAllele=unlist(alt, use.names=FALSE), ...) 
        origID <- rep(seq_len(length(rd)), elementLengths(alt))
        res$queryID <- origID[res$queryID]
        res
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
        on.exit(isActiveSeq(subject) <- masks)
        .setActiveSubjectSeq(query, subject)

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

        ## map genome position to transcript position
        cdsByTx <- cdsByTx[unique(subjectHits(fo))]
        txlocal <- globalToLocal(queryAdj, cdsByTx)
        midx <- which(width(varAllele) == 0)
        if (length(midx) > 0) {
            warning("records with missing 'varAllele' values ",
                    "will be ignored")
            txlocal <- txlocal[!txlocal$qindex %in% midx, ] 
        }
        xCoding <- query[txlocal$qindex]
        xAllele <- varAllele[txlocal$qindex]
 
        ## construct original sequences 
        originalWidth <- width(xCoding)
        codonStart <- ((start(txlocal$txloc) - 1L) %/% 3L) * 3L + 1L
        ## codonEnd must be adjusted for 
        ## (1) the width of the reference sequence and
        ## (2) the position in the codon of the alternate allele substitution
        varPosition <- (start(txlocal$txloc) - 1L) %% 3L + 1L
        codonEnd <- 
            codonStart + (((varPosition + originalWidth) %/% 3L) * 3L + 2L)
        txseqs <- getTranscriptSeqs(cdsByTx, seqSource)
        codons <- DNAStringSet(substring(txseqs[txlocal$sindex], 
            codonStart, codonEnd))

        ## construct variant sequences 
        varWidth <- width(xAllele)
        indels <- originalWidth == 0 | varWidth == 0
        translateIdx <- abs(varWidth - originalWidth) %% 3 == 0 
        n <- grep("N", as.character(xAllele, use.names=FALSE), fixed=TRUE)
        if (length(n) > 0) {
            warning("varAllele values containing 'N' will not be translated")
            translateIdx[n] <- FALSE
        }
        varSeq <- codons
        subseq(varSeq, start=varPosition, width=originalWidth) <- xAllele 

        ## results
        queryID <- txlocal$qindex
        txLoc <- txlocal$txloc
        txID <- names(cdsByTx)[txlocal$sindex]
        cdsID <- txlocal$cdsid 
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
 
        DataFrame(queryID, txLoc, consequence, refSeq, varSeq, 
            refAA, varAA, txID, cdsID, geneID) 
    }
)

