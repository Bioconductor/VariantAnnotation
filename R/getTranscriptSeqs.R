getTranscriptSeqs <- function(cds_by_tx, BSgenomeOrganism) 
{
    #cds_seqs <- getSeq(organism, unlist(cds_by_tx, use.names=FALSE))
    #tx_id <- rep(seq(length(cds_by_tx)), elementLengths(cds_by_tx))
    #unlist(lapply(split(cds_seqs, tx_id), paste, collapse=""))

    ## BSgenome 
    extractTranscriptsFromGenome(BSgenomeOrganism, cds_by_tx)

    ## bam file approach
    ## overhead :
    ## - need repository of files, iterate through files for each chrom
    ## - want one seq per tx composed of cds only
    ##   need to id ranges of all cds per tx then combine 
    #fl <- "~/proj/data/chr16.subst.fa"
    ##indexFa(fl) 
    #fa <- open(FaFile(fl))
    #countFa(fa)
    ##idx <- scanFaIndex(fa)
    #idx <- cds_by_tx[40000] 
    #dna <- scanFa(fa, unlist(idx))
    ##ranges(idx) <- narrow(range(idx), -10)
}

