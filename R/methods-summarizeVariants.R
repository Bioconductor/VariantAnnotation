### =========================================================================
### summarizeVariants methods 
### =========================================================================

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "CodingVariants"),
    function(subject, query, mode, ...)
{
    grl <- cdsBy(subject, "tx")
    callGeneric(grl, query, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "FiveUTRVariants"),
    function(subject, query, mode, ...)
{
    grl <- fiveUTRsByTranscript(subject) 
    callGeneric(grl, query, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "ThreeUTRVariants"),
    function(subject, query, mode, ...)
{
    grl <- threeUTRsByTranscript(subject) 
    callGeneric(grl, query, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "SpliceSiteVariants"),
    function(subject, query, mode, ...)
{
    grl <- intronsByTranscript(subject) 
    callGeneric(grl, query, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "IntronVariants"),
    function(subject, query, mode, ...)
{
    grl <- intronsByTranscript(subject) 
    callGeneric(grl, query, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "PromoterVariants"),
    function(subject, query, mode, ...)
{
    gr <- transcripts(subject, columns="tx_id")
    grl <- splitAsList(gr, seq_len(length(gr))) 
    names(grl) <- mcols(gr)$tx_id
    callGeneric(grl, query, mode, ...)
})

setMethod("summarizeVariants", c("GRangesList", "VCF", "ANY"),
    function(subject, query, mode, ..., subjectFactor=factor(seq_len(length(subject))))
{
    if (length(subjectFactor) != length(subject))
        stop("'subjectFactor' must be the same length as 'subject'")
    if (any(na <- is.na(subjectFactor))) {
        warning("NA 'subjectFactor' levels were removed")
        subject <- subject[!na]
        subjectFactor <- subjectFactor[!na]
    } 
    if (!is.factor(subjectFactor))
        subjectFactor <- factor(subjectFactor)
    if (length(geno(query)) == 0L) {
        warning("No genotypes found in 'query'.")
        return(.baseSE(query, subject))
    }

    if (is.function(mode)) {
    ## findOverlaps() 
        ct <- mode(rowData(query), subject, ...)
        if (length(ct) == 0L)
            return(.baseSE(query, subject))
        hits <- unique(data.frame(queryHits(ct),
                       subjectFactor[subjectHits(ct)])) 
    } else {
    ## locateVariants() 
        ct <- locateVariants(rowData(query), subject, mode, ...)
        if (length(ct) == 0L)
            return(.baseSE(query, subject))
        subid <- match(ct$TXID, names(subject))
        hits <- unique(data.frame(ct$QUERYID, subjectFactor[subid]))
    }

    ## relist data to reflect subjectFactor
    if (identical(subjectFactor, factor(seq_len(length(subject)))))
        rd <- subject 
    else
        rd <- splitAsList(subject@unlistData, 
                          rep(subjectFactor, elementLengths(subject)))

    ## genotypes
    na <- c("0|0", "0/0", "./.", ".|.")
    gtype <- as.numeric(!(geno(query)$GT[unique(hits[,1]),]) %in% na)

    ## summarize variants factor-by-sample 
    facbyvar <- table(hits[,2], hits[,1])
    varbysamp <- matrix(gtype, ncol=ncol(query))
    facbysamp <- facbyvar %*% varbysamp

    SummarizedExperiment(rowData=rd, colData=colData(query), 
                         exptData=exptData(query),
                         assays=SimpleList(counts=facbysamp))
})

.baseSE <- function(query, subject, ...)
{
    SummarizedExperiment(
        rowData=subject,
        colData=colData(query),
        exptData=exptData(query),
        assays=SimpleList(counts=matrix(NA_integer_, 
            nrow=length(subject), ncol=ncol(query))))
}
