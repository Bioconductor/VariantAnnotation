annotateVariants <- 
    function(query, subject, orgPackage, BSgenomeOrganism, variantColumn, ...)
{

    ## location
    loc <- locateVariants(query, subject, orgPackage) 
 
    ## protein coding changes
    proteins <- predictCoding(query, subject, BSgenomeOrganism, variantColumn) 
 
    aa <- substr(values(proteins)[["aa"]], start=1, stop=length(proteins))
    aa.ref <- substr(values(proteins)[["aa.ref"]], start=1, stop=length(proteins))
    co <- countOverlaps(query, proteins)
    queryHit <- co[co != 0]
    AAcoding <- vector("list", length(query))
    AAcoding[co != 0] <- as.list(split(values(proteins)[c("tx_id", "aa", "aa.ref")], 
        rep(seq_len(length(queryHit)), queryHit)))

    loc$AAcoding <- AAcoding 
    values(query) <- append(values(query), loc)
    query
}

