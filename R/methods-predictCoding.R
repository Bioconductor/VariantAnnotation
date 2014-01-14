### =========================================================================
### predictCoding methods 
### =========================================================================

setMethod("predictCoding", c("Ranges", "TranscriptDb", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    callGeneric(as(query, "GRanges"), subject, seqSource, varAllele, ...,
                ignore.strand=ignore.strand) 
})

setMethod("predictCoding", c("CollapsedVCF", "TranscriptDb", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    rd <- rowData(query)
    alt <- alt(query) 
    if (is(alt, "CharacterList")) {
        alt <- .toDNAStringSetList(alt)
        if (sum(elementLengths(alt)) == 0L) {
            stop("No nucleotide ALT values were detected.")
        }
    }
    rd <- rep(rowData(query), elementLengths(alt))
    res <- callGeneric(rd, subject, seqSource, unlist(alt, use.names=FALSE), 
                ..., ignore.strand=ignore.strand)
    ## adjust QUERYID for expansion of rowData
    res$QUERYID <- rep(seq_len(length(alt)),
                       elementLengths(alt))[res$QUERYID]
    res 
})

setMethod("predictCoding", c("ExpandedVCF", "TranscriptDb", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    ## ExpandedVCF ALT must be DNAStringSet (CharacterList not supported)
    callGeneric(rowData(query), subject, seqSource, alt, ..., 
                ignore.strand=ignore.strand) 
})

setMethod("predictCoding", c("GRanges", "TranscriptDb", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    .predictCoding(query, subject, seqSource, varAllele, ...,
                   ignore.strand=ignore.strand)
})

.predictCoding <-
    function(query, subject, seqSource, varAllele, ..., 
             cache=new.env(parent=emptyenv()), ignore.strand=FALSE)
{
    stopifnot(length(varAllele) == length(query))
    if (!any(seqlevels(query) %in% seqlevels(subject)))
        warning("none of seqlevels(query) match seqlevels(subject)")

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

    txlocal <- .predictCodingGRangesList(query, cache[["cdsbytx"]], seqSource, 
                                         varAllele, ignore.strand=ignore.strand)
    txid <- values(txlocal)[["TXID"]] 
    values(txlocal)[["GENEID"]] <- map$geneid[match(txid, map$txid)]
    txlocal
}

.predictCodingGRangesList <- function(query, cdsbytx, seqSource, varAllele, ...,
                                      ignore.strand=FALSE)
{
    ## FIXME : set query back after olaps
    if (any(insertion <- width(query) == 0))
        start(query)[insertion] <- start(query)[insertion] - 1
    if (ignore.strand)
        strand(query) <- "*"

    ## retrieve local coordinates
    values(query) <- append(values(query), DataFrame(varAllele=varAllele))
    txlocal <- refLocsToLocalLocs(ranges=query, cdsbytx=cdsbytx)

    if (length(txlocal) == 0)
        return(txlocal)
 
    rwidth <- width(txlocal)
    translateidx <- rep(TRUE, length(txlocal)) 
    ## reverse complement variant alleles for "-" strand
    nstrand <- as.vector(strand(txlocal) == "-")
    if (any(nstrand))
        values(txlocal)[["varAllele"]][nstrand] <-
            reverseComplement(values(txlocal)[["varAllele"]][translateidx & nstrand])
    altallele <- values(txlocal)[["varAllele"]]
    fmshift <- abs(width(altallele) - rwidth) %% 3 != 0 
    if (any(fmshift))
        translateidx[fmshift] <- FALSE
    zwidth <- width(altallele) == 0
    if (any(zwidth)) {
        warning("records with missing 'varAllele' were ignored")
        translateidx[zwidth] <- FALSE 
        fmshift[zwidth] <- FALSE
    }
    codeN <- rep(FALSE, length(txlocal)) 
    codeN[grep("N", as.character(altallele, use.names=FALSE), 
        fixed=TRUE)] <-TRUE
    if (any(codeN)) {
        warning("varAllele values containing 'N' were not translated")
        translateidx[codeN] <- FALSE
    }

    ## reference and variant codon sequences
    altpos <- (start(values(txlocal)[["CDSLOC"]]) - 1L) %% 3L + 1L
    refCodon <- varCodon <- .constructRefSequences(txlocal, altpos, seqSource, 
                                                   cdsbytx)
    if (any(translateidx)) {
        subseq(varCodon, altpos, width=rwidth)[translateidx]  <-
            altallele[translateidx]

        ## translation
        refAA <- translate(refCodon)
        varAA <- AAStringSet(rep("", length(txlocal))) 
        varAA[translateidx] <- translate(varCodon[translateidx])
    } else {
        refAA <- varAA <- AAStringSet(rep("", length(txlocal))) 
    }

    ## results
    nonsynonymous <- as.character(refAA) != as.character(varAA) 
    consequence <- rep("synonymous", length(txlocal))
    consequence[nonsynonymous] <- "nonsynonymous" 
    consequence[fmshift] <- "frameshift"
    consequence[nonsynonymous & (as.character(varAA) %in% "*")] <- "nonsense" 
    consequence[zwidth | codeN] <- "not translated" 
    consequence <- factor(consequence) 
 
    values(txlocal) <- append(values(txlocal), 
        DataFrame(GENEID=NA_character_, 
                  CONSEQUENCE=consequence, 
                  REFCODON=refCodon, 
                  VARCODON=varCodon, 
                  REFAA=refAA, VARAA=varAA))
    txlocal 
}

.constructRefSequences <- function(txlocal, altpos, seqSource, cdsbytx)
{ 
    ## adjust codon end for 
    ## - width of the reference sequence
    ## - position of alt allele substitution in the codon
    cstart <- ((start(values(txlocal)[["CDSLOC"]]) - 1L) %/% 3L) * 3L + 1L
    cend <- cstart + (((altpos + width(txlocal) - 2L) %/% 3L) * 3L + 2L)
    txord <- match(values(txlocal)[["TXID"]], names(cdsbytx))
    txseqs <- getTranscriptSeqs(cdsbytx[txord], seqSource)
    DNAStringSet(substring(txseqs, cstart, cend))
}
