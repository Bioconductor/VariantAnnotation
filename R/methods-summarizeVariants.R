### =========================================================================
### summarizeVariants methods 
### =========================================================================

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "CodingVariants"),
    function(query, subject, mode, ...)
{
    grl <- cdsBy(query, "tx")
    callGeneric(grl, subject, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "FiveUTRVariants"),
    function(query, subject, mode, ...)
{
    grl <- fiveUTRsByTranscript(query) 
    callGeneric(grl, subject, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "ThreeUTRVariants"),
    function(query, subject, mode, ...)
{
    grl <- threeUTRsByTranscript(query) 
    callGeneric(grl, subject, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "SpliceSiteVariants"),
    function(query, subject, mode, ...)
{
    grl <- intronsByTranscript(query) 
    callGeneric(grl, subject, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "IntronVariants"),
    function(query, subject, mode, ...)
{
    grl <- intronsByTranscript(query) 
    callGeneric(grl, subject, mode, ...)
})

setMethod("summarizeVariants", c("TranscriptDb", "VCF", "PromoterVariants"),
    function(query, subject, mode, ...)
{
    gr <- transcripts(query, columns="tx_id")
    grl <- splitAsList(gr, seq_len(length(gr))) 
    names(grl) <- mcols(gr)$tx_id
    callGeneric(grl, subject, mode, ...)
})

setMethod("summarizeVariants", c("GRangesList", "VCF", "VariantType"),
    function(query, subject, mode, ...)
{
    callGeneric(query, subject, mode=locateVariants, ..., region=mode, 
        asHits=TRUE)
})

setMethod("summarizeVariants", c("GRangesList", "VCF", "function"),
    function(query, subject, mode, ...)
{
    if (length(geno(subject)) == 0L) {
        warning("No genotypes found in 'query'.")
        return(.baseSE(query, subject))
    }
    ## count
    hits <- mode(rowData(subject), query, ...)
        if (length(hits) == 0L)
            return(.baseSE(query, subject))

    ## genotypes
    na <- c("0|0", "0/0", "./.", ".|.")
    vcf_geno <- geno(subject)$GT[unique(queryHits(hits)), ]
    gtype <- as.numeric(!vcf_geno %in% na)

    ## summarize counts factor-by-sample 
    fac_x_var <- table(subjectHits(hits), queryHits(hits))
    var_x_smp <- matrix(gtype, ncol=ncol(subject))
    fac_x_smp <- fac_x_var %*% var_x_smp

    SummarizedExperiment(rowData=query[unique(subjectHits(hits))], 
                         colData=colData(subject), 
                         exptData=exptData(subject),
                         assays=SimpleList(counts=fac_x_smp))
})

.baseSE <- function(query, subject, ...)
{
    SummarizedExperiment(rowData=query, colData=colData(subject),
                         exptData=exptData(subject),
                         assays=SimpleList(counts=matrix(NA_integer_, 
                             nrow=length(query), ncol=ncol(subject))))
}
