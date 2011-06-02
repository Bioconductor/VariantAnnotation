getTranscriptSeqs <- function(cds_by_tx, BSgenomeOrganism) 
{
    extractTranscriptsFromGenome(BSgenomeOrganism, cds_by_tx)
}

