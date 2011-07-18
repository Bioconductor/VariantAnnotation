## FIXME : (1) overlaps all type = 'within' or 'any' ?
##         (2) if doesn't match a tx -> classified as intergenic 
 
setMethod("locateVariants", c("GRanges", "TranscriptDb"),
    function(query, subject, ...)
    {
        chrom <- seqlevels(query)
        isActiveSeq(subject)[!names(isActiveSeq(subject)) %in% chrom] <- FALSE 
        cdsByTx <- cdsBy(subject)
        tx <- transcripts(subject, columns=c("exon_id", "tx_id", "gene_id"))

        ## adjust query with width=0
        ## - increment end to equal start value 
        ## - findOverlaps - findOverlaps(type="start") 
        if (any(width(query) == 0)) {
            insertion <- width(query) == 0
            queryAdj <- query
            start(queryAdj[insertion]) <- 
                start(query)[insertion] - 1
        } else queryAdj <- query

### insertions in coding regions
#foany <- findOverlaps(queryAdj, cdsByTx)
#fostart <- findOverlaps(queryAdj, cdsByTx, type="start")
#m1 <- matchMatrix(foany)
#m2 <- matchMatrix(fostart)
#key <- function(m)paste(sep="\1",m[,1],m[,2])
#m1[!is.element(key(m1),key(m2)),]
#
### insertions in transcripts 
#foany <- findOverlaps(queryAdj, tx)
#fostart <- findOverlaps(queryAdj, tx, type="start")
#m1 <- matchMatrix(foany)
#m2 <- matchMatrix(fostart)
#key <- function(m)paste(sep="\1",m[,1],m[,2])
#m1[!is.element(key(m1),key(m2)),]


        cdsFO <- findOverlaps(queryAdj, cdsByTx, type="within")
        txFO <- findOverlaps(queryAdj, tx, type="within")
        cdsCO <- tabulate(queryHits(cdsFO), length(query))
        txCO <- tabulate(queryHits(txFO), length(query))

        queryIndex <- queryHits(txFO)
        txID <- names(cdsByTx)[subjectHits(txFO)]
        geneID <- values(tx)[["gene_id"]][subjectHits(txFO)]

        ## coding :
        coding <- cdsCO > 0 

        ## intron :
        intron <- txCO != 0 & cdsCO == 0

t1 <- proc.time() 
#splicings <- GenomicFeatures:::.getFullSplicings(subject,
#    translated.transcripts.only=TRUE)
#fiveUTR <- my5UTR(subject, splicings)
#threeUTR <- my3UTR(subject, splicings)

        ## UTRs :
        fiveUTR <- fiveUTRsByTranscript(subject)
        threeUTR <- threeUTRsByTranscript(subject)
t2 <- proc.time() 

cat("diff :", t2 - t1, "\n")

        utr5 <- countOverlaps(queryAdj, fiveUTR, type="within") > 0
        utr3 <- countOverlaps(queryAdj, threeUTR, type="within") > 0

        Location <- rep("unknown", length(queryIndex))
        Location[queryIndex %in% which(intron)] <- "intron"
        Location[queryIndex %in% which(utr5)] <- "5'UTR"
        Location[queryIndex %in% which(utr3)] <- "3'UTR"
        Location[queryIndex %in% which(coding)] <- "coding"

        dat1 <- DataFrame(queryHits=queryIndex, txID, geneID, Location)

        ## intergenic :
        ## FIXME : ridiculous, need to simplify
        if (any(txCO == 0)) {
            intergenic <- txCO == 0
            intvar <- query[txCO == 0]
            flankGenes <- lapply(seq_len(length(intvar)), 
                function(x, tx, intvar) {
                    pre <- fol <- intvar[x,]
                    pregene <- values(tx[precede(pre, tx)])["gene_id"] 
                    while (length(pregene$gene_id[[1]]) == 0 && 
                           (typeof(pregene$gene_id[[1]]) == "character")){
                            pre <- precede(pre, tx)
                            pregene <- values(tx[precede(pre, tx)])["gene_id"] 
                    }

                    folgene <- values(tx[follow(fol, tx)])["gene_id"] 
                    while (length(folgene$gene_id[[1]]) == 0 && 
                           (typeof(folgene$gene_id[[1]]) == "character")){
                            fol <- tx[follow(fol, tx)]
                            folgene <- values(tx[follow(fol, tx)])["gene_id"] 
                    }
                c(pregene$gene_id[[1]], folgene$gene_id[[1]])
                }, tx, intvar)
        } else flankGenes <- NA

       dat2 <- DataFrame(queryHits=which(intergenic), 
           txID=rep(NA, length(which(intergenic))), 
           geneID=CharacterList(flankGenes), 
           Location=rep("intergenic", length(which(intergenic))))

       ans <- rbind(dat1, dat2)
       ans[order(ans$queryIndex), ]
    }
)


