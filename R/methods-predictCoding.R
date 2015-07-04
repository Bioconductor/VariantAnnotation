### =========================================================================
### predictCoding methods 
### =========================================================================

setMethod("predictCoding", c("Ranges", "TxDb", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    callGeneric(as(query, "GRanges"), subject, seqSource, varAllele, ...,
                ignore.strand=ignore.strand) 
})

setMethod("predictCoding", c("CollapsedVCF", "TxDb", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    rd <- rowRanges(query)
    alt <- alt(query) 
    if (is(alt, "CharacterList")) {
        alt <- .toDNAStringSetList(alt)
        if (sum(elementLengths(alt)) == 0L) {
            stop("No nucleotide ALT values were detected.")
        }
    }
    rd <- rep(rowRanges(query), elementLengths(alt))
    res <- callGeneric(rd, subject, seqSource, unlist(alt, use.names=FALSE), 
                ..., ignore.strand=ignore.strand)
    ## adjust QUERYID for expansion of rowRanges
    res$QUERYID <- rep(seq_len(length(alt)),
                       elementLengths(alt))[res$QUERYID]
    res 
})

setMethod("predictCoding", c("ExpandedVCF", "TxDb", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    if (is(alt(query), "CharacterList")) {
      stop("alt(query) must be a DNAStringSet (not a CharacterList)")
    }
    callGeneric(rowRanges(query), subject, seqSource, alt(query), ..., 
                ignore.strand=ignore.strand) 
})

setMethod("predictCoding", c("GRanges", "TxDb", "ANY", "DNAStringSet"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    .predictCoding(query, subject, seqSource, varAllele, ...,
                   ignore.strand=ignore.strand)
})

setMethod("predictCoding", c("VRanges", "TxDb", "ANY", "missing"),
    function(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)
{
    varAllele <- alt(query)
    query <- as(query, "GRanges")
    if (!is(varAllele, "DNAStringSet")) {
        tryCatch({
            varAllele <- DNAStringSet(varAllele)
        }, error=function(e) {
            stop(paste0("attempt to coerce 'alt' to DNAStringSet failed with ",
                 "error: ", conditionMessage(e)))
        })
    }
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
                      txid=mcols(unlist(cache[["txbygene"]], 
                          use.names=FALSE))[["tx_id"]],
                      stringsAsFactors=FALSE)

    txlocal <- .predictCodingGRangesList(query, cache[["cdsbytx"]], 
                   seqSource, varAllele, ..., ignore.strand=ignore.strand)
    txid <- mcols(txlocal)$TXID 
    mcols(txlocal)$GENEID <- map$geneid[match(txid, map$txid)]
    txlocal
}

.predictCodingGRangesList <- function(query, cdsbytx, seqSource, varAllele, 
                                      ..., genetic.code=GENETIC_CODE,
                                      if.fuzzy.codon="error", 
                                      ignore.strand=FALSE)
{
    if (ignore.strand)
        strand(query) <- "*"

    ## variant location in cds region
    mcols(query) <- append(mcols(query), DataFrame(varAllele=varAllele))
    txlocal <- .localCoordinates(query, cdsbytx, ignore.strand=FALSE, ...)
    if (length(txlocal) == 0)
        return(txlocal)

    ## reverse complement "-" strand
    valid <- rep(TRUE, length(txlocal))
    nstrand <- as.vector(strand(txlocal) == "-")
    if (any(nstrand)) {
        va <- mcols(txlocal)$varAllele
        va[nstrand] <- reverseComplement(va[valid & nstrand])
        mcols(txlocal)$varAllele <- va
    }

    ## frameshift
    refwidth <- width(txlocal)
    altallele <- mcols(txlocal)$varAllele
    fmshift <- abs(width(altallele) - refwidth) %% 3 != 0 
    if (any(fmshift))
        valid[fmshift] <- FALSE

    ## zero-width
    zwidth <- width(altallele) == 0
    if (any(zwidth)) {
        warning("records with missing 'varAllele' were ignored")
        valid[zwidth] <- FALSE 
        fmshift[zwidth] <- FALSE
    }

    ## reference codon sequences
    altpos <- (start(mcols(txlocal)$CDSLOC) - 1L) %% 3L + 1L
    refCodon <- varCodon <- .getRefCodons(txlocal, altpos, seqSource, cdsbytx)

    ## allowed characters that can't be translated
    ## "N", ".", "+" and "-"
    pattern <- "N|\\.|\\+|\\-"
    altCheck <- grepl(pattern, as.character(altallele, use.names=FALSE))
    refCheck <- grepl(pattern, as.character(refCodon, use.names=FALSE))
    noTrans <- rep(FALSE, length(txlocal)) 
    noTrans[altCheck | refCheck] <- TRUE
    valid[noTrans] <- FALSE
    if (any(altCheck))
        warning("'varAllele' values with 'N', '.', '+' or '-'",
                " were not translated")
    if (any(refCheck))
        warning("reference codons with 'N', '.', '+' or '-'",
                " were not translated")

    ## substitute and translate
    refAA <- varAA <- AAStringSet(rep("", length(txlocal))) 
    if (any(valid)) {
        subseq(varCodon, altpos, width=refwidth)[valid] <- altallele[valid]
        refAA[valid] <- translate(refCodon[valid], genetic.code=genetic.code,
                                  if.fuzzy.codon=if.fuzzy.codon)
        varAA <- AAStringSet(rep("", length(txlocal))) 
        varAA[valid] <- translate(varCodon[valid], genetic.code=genetic.code, 
                                  if.fuzzy.codon=if.fuzzy.codon)
    }
    varCodon[!valid] <- ""

    ## results
    nonsynonymous <- as.character(refAA) != as.character(varAA) 
    consequence <- rep("synonymous", length(txlocal))
    consequence[nonsynonymous] <- "nonsynonymous" 
    consequence[fmshift] <- "frameshift"
    consequence[nonsynonymous & (as.character(varAA) %in% "*")] <- "nonsense" 
    consequence[zwidth | noTrans] <- "not translated" 
    consequence <- factor(consequence) 
 
    mcols(txlocal) <- append(mcols(txlocal), 
        DataFrame(GENEID=NA_character_, 
                  CONSEQUENCE=consequence, 
                  REFCODON=refCodon, 
                  VARCODON=varCodon, 
                  REFAA=refAA, VARAA=varAA))
    txlocal 
}

.getRefCodons <- function(txlocal, altpos, seqSource, cdsbytx)
{ 
    ## adjust codon end for 
    ## - width of the reference sequence
    ## - position of alt allele substitution in the codon
    cstart <- ((start(mcols(txlocal)$CDSLOC) - 1L) %/% 3L) * 3L + 1L
    cend <- cstart + (((altpos + width(txlocal) - 2L) %/% 3L) * 3L + 2L)
    txord <- match(mcols(txlocal)$TXID, names(cdsbytx))
    txseqs <- extractTranscriptSeqs(seqSource, cdsbytx[txord])
    DNAStringSet(substring(txseqs, cstart, cend))
}
