
setMethod("getTranscriptSeqs",  c("BSgenome", "GRangesList"),
    function(seqSource, subject, ...)
    {
    extractTranscriptsFromGenome(seqSource, subject)
    }
)

setMethod("getTranscriptSeqs",  c("FaFileList", "GRangesList"),
    function(seqSource, subject, ...)
    {
    faseqs <- lapply(seqSource, function(x) {
        ## FIXME : close x each time? 
        seqlevels(scanFaIndex(open(x)))
        })
    seqmatch <- seqlevels(subject) %in% faseqs 
    if (any(seqmatch == FALSE)){
        missingSeq <- seqlevels(subject)[seqmatch == FALSE]
        stop(paste("sequence ", missingSeq, 
            " not found in 'seqSource' \n", sep=""))
    }
    usubj <- unlist(subject, use.names=FALSE) 
    txmap <- rep(names(subject), elementLengths(subject))
    chrmap <- rep(runValue(seqnames(usubj)),
        runLength(seqnames(usubj)))
    resultmap <- rep(seq_len(length(subject)), elementLengths(subject))
    widths <- sum(width(subject)) 

    ## FIXME : handle strand
    lst <- sapply(seqSource, function(fl, usubj, chrmap, txmap) {
        fa <- open(fl)
        idx <- scanFaIndex(fa)
        matchidx <- chrmap == seqlevels(idx)
        rng <- usubj[matchidx] 
        dna <- scanFa(fa, rng)
        close(fa)
 
        vwidth <- widths[unique(resultmap[matchidx])]
        udna <- unlist(dna)
        vws <- successiveViews(udna, vwidth)
        tx <- as(vws, "DNAStringSet") 
        names(tx) <- unique(txmap[matchidx])
        metadata(tx) <- list(index=unique(resultmap[matchidx])) 
        tx
    }, usubj, chrmap, txmap, USE.NAMES=FALSE)

    ## reorder
    ans <- rep.int(as("", "DNAStringSet"), length(subject)) 
    all <- unlist.list.of.XVectorList("DNAStringSet", lst)
    orderIdx <- unlist(lapply(lst, function(x) metadata(x)), use.names=FALSE)
    ans <- all[orderIdx]
    ans
    }
)

