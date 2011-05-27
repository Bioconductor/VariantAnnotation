#predictCoding <- function(x, cds_by_tx, organism=Hsapiens, 
#    tx_seqs = getTranscriptSeqs(cds_by_tx, organism), variant_column = "read")
predictCoding <- 
    function(x, txdb, BSgenomeOrganism = Hsapiens, variantColumn = "read")
{
    cds_by_tx <- cdsBy(txdb)
    tx_seqs <- getTranscriptSeqs(cds_by_tx, BSgenomeOrganism) 
    x <- as(x, "GRanges")
    tx_local <- globalToLocal(x, cds_by_tx)
    x_coding <- x[tx_local$global.ind]

    ## extract codon sequences
    codon_start <- (start(tx_local$local) - 1L) %/% 3L * 3L + 1L
    codons <- substring(tx_seqs[tx_local$ranges.ind], codon_start, 
        codon_start + 2L)
    codon_pos <- (start(tx_local$local) - 1L) %% 3L + 1L

    ## assuming variant with single base pair change (ie, snv)
    var_codons <- codons
    substring(var_codons, codon_pos, codon_pos) <- 
        as.character(values(x_coding)[[variantColumn]])

    values(x_coding)[colnames(values(cds_by_tx))] <-
        values(cds_by_tx)[tx_local$ranges.ind,]
    values(x_coding)$tx_id <- names(cds_by_tx)[tx_local$ranges.ind]
    values(x_coding)$aa <- translate(DNAStringSet(var_codons))
    values(x_coding)$aa.ref <- translate(DNAStringSet(codons))
    x_coding
}

