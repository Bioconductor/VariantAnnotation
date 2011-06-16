
setMethod("getTranscriptSeqs",  c("GRangesList", "BSgenome"),
    function(query, subject, ...)
    {
    extractTranscriptsFromGenome(subject, query)
    }
)

setMethod("getTranscriptSeqs",  c("GRangesList", "FaFile"),
    function(query, subject, ...)
    {
        uquery <- unlist(query, use.names=FALSE)
        widths <- sum(width(query)) 
        seq <- callGeneric(uquery, subject, ...)
        vws <- successiveViews(seq, widths)
        tx <- as(vws, "DNAStringSet") 
        names(tx) <- names(query) 
        tx
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
        dna <- scanFa(fa, query)
        ans <- unlist(dna)
        close(fa)
        ans 
 
        ## FIXME : handle strand
    }
)

