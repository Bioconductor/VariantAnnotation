### =========================================================================
### getTranscriptSeqs methods 
### =========================================================================

setMethod("getTranscriptSeqs",  c("GRangesList", "BSgenome"),
    function(query, subject, ...)
    {
    extractTranscriptsFromGenome(subject, query)
    }
)

setMethod("getTranscriptSeqs",  c("GRangesList", "FaFile"),
    function(query, subject, ...)
    {
        ## order by exon rank, check for mixed chromosomes, mixed strand 
        txlist <- GenomicFeatures:::.makeUCSCTxListFromGRangesList(query,
            decreasing.rank.on.minus.strand=FALSE)

        ## back to GRanges for getSeq method 
        widthEL <- elementLengths(query)
        gr <- GRanges(
                seqnames=Rle(txlist$chrom, widthEL), 
                ranges=IRanges(start=unlist(txlist$exonStarts), 
                               end=unlist(txlist$exonEnds)), 
                strand=Rle(txlist$strand, widthEL))

        seqs <- callGeneric(gr, subject, ...)
        oneseq <- unlist(seqs)
        widthView <- sum(width(query)) 
        vws <- successiveViews(oneseq, widthView)
        txseq <- as(vws, "DNAStringSet") 
        names(txseq) <- names(query) 
        txseq
    }
)
 
setMethod("getTranscriptSeqs",  c("GRanges", "FaFile"),
    function(query, subject, ...)
    {
        fa <- open(subject)
        idx <- scanFaIndex(fa)
        seqmatch <- seqlevels(query) %in% seqlevels(idx) 
        if (any(seqmatch == FALSE)){
            missingSeq <- seqlevels(subject)[seqmatch == FALSE]
            stop(paste("sequence ", missingSeq, 
                " not found in 'query' \n", sep=""))
        }
 
        dna <- getSeq(fa, query) 
        close(fa)
        dna 
    }
)

