### =========================================================================
### predictCoding methods 
### =========================================================================

setMethod("predictCoding",  
    signature("Ranges", "TranscriptDb", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ...)
    {
        x <- as(query, "GRanges")
        callGeneric(query=x, subject=subject, seqSource, 
            varAllele, ...) 
    }
)

setMethod("predictCoding", 
    signature("VCF", "TranscriptDb", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ...)
    {
        alt <- values(alt(query))[["ALT"]]
        if (!is(alt, "DNAStringSetList"))
            stop("alt(<VCF>)[['ALT']] must be a DNAStringSetList")
        rd <- rep(rowData(query), elementLengths(alt))
        callGeneric(query=rd, subject=subject, seqSource=seqSource, 
            varAllele=unlist(alt, use.names=FALSE), ...) 
    }
)

setMethod("predictCoding", 
    signature("GRanges", "TranscriptDb", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ...)
{
    .predictCoding(query, subject, seqSource, varAllele, ...)
})

.predictCoding <-
    function(query, subject, seqSource, varAllele, ..., 
             cache=new.env(parent=emptyenv()))
{
    stopifnot(length(varAllele) == length(query))
    if (!any(seqlevels(query) %in% seqlevels(subject)))
        warning("none of seqlevels(query) match seqlevels(subject)")

    ## mask chromosomes not in query
    masks <- isActiveSeq(subject)
    on.exit(isActiveSeq(subject) <- masks)
    .setActiveSubjectSeq(query, subject)

    if (!exists(".__init__", cache, inherits=FALSE)) {
        cache[["cdsbytx"]] <- cdsBy(subject)
        cache[["txbygene"]] <- transcriptsBy(subject, "gene")
        cache[[".__init__"]] <- TRUE
    }

    map <- data.frame(geneid=rep(names(cache[["txbygene"]]), 
                          elementLengths(cache[["txbygene"]])),
                      txid=values(unlist(cache[["txbygene"]], 
                          use.names=FALSE))[["tx_id"]],
                      stringsAsFactors=FALSE)

    ## FIXME : set query back after olaps
    if (any(insertion <- width(query) == 0))
        start(query)[insertion] <- start(query)[insertion] - 1

    ## retrieve local coordinates
    values(query) <- append(values(query), DataFrame(varAllele=varAllele))
    txlocal <- refLocsToLocalLocs(ranges=query, cdsbytx=cache[["cdsbytx"]])

    if (length(txlocal) == 0)
        return(txlocal)
    rwidth <- width(txlocal)
    translateidx <- rep(TRUE, length(txlocal)) 
    altallele <- values(txlocal)[["varAllele"]]
    fmshift <- abs(width(altallele) - rwidth) %% 3 != 0 
    if (any(fmshift))
        translateidx[fmshift] <- FALSE
    zwidth <- width(altallele) == 0
    if (any(zwidth)) {
        warning("records with missing 'varAllele' values will be ignored")
        translateidx[zwidth] <- FALSE 
        fmshift[zwidth] <- FALSE
    }
    codeN <- rep(FALSE, length(txlocal)) 
    codeN[grep("N", as.character(altallele, use.names=FALSE), 
        fixed=TRUE)] <-TRUE
    if (any(codeN)) {
        warning("varAllele values containing 'N' will not be translated")
        translateidx[codeN] <- FALSE
    }

    ## reference and variant codon sequences
    altpos <- (start(values(txlocal)[["cdsLoc"]]) - 1L) %% 3L + 1L
    refCodon <- varCodon <- .constructRefSequences(txlocal, altpos, seqSource, cache)
    subseq(varCodon, start=altpos, width=rwidth) <- altallele

    ## translation
    refAA <- translate(refCodon)
    varAA <- AAStringSet(rep("", length(txlocal))) 
    varAA[translateidx] <- translate(varCodon[translateidx])

    ## results
    txID <- values(txlocal)[["txID"]] 
    geneID <- map$geneid[match(txID, map$txid)]
    nonsynonymous <- as.character(refAA) != as.character(varAA) 
    consequence <- rep("synonymous", length(txlocal))
    consequence[nonsynonymous] <- "nonsynonymous" 
    consequence[fmshift] <- "frameshift" 
    consequence[zwidth | codeN] <- "not translated" 
    consequence <- factor(consequence) 
 
    values(txlocal) <- append(values(txlocal), DataFrame(geneID, consequence, 
        refCodon, varCodon, refAA, varAA))
    txlocal 
}

.constructRefSequences <- function(txlocal, altpos, seqSource, cache)
{ 
    ## adjust codon end for 
    ## - width of the reference sequence
    ## - position of alt allele substitution in the codon
    cstart <- ((start(values(txlocal)[["cdsLoc"]]) - 1L) %/% 3L) * 3L + 1L
    cend <- cstart + (((altpos + width(txlocal) - 2L) %/% 3L) * 3L + 2L)
    txord <- match(values(txlocal)[["txID"]], names(cache[["cdsbytx"]]))
    txseqs <- getTranscriptSeqs(cache[["cdsbytx"]][txord], seqSource)
    DNAStringSet(substring(txseqs, cstart, cend))
}
