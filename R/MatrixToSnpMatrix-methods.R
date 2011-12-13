### =========================================================================
### MatrixToSnpMatrix methods 
### =========================================================================

setMethod("MatrixToSnpMatrix", c("matrix", "GRanges"),
    function(callMatrix, mapGRanges, ...)
    {
        ## Coding for snpMatrix :
        ## 0 = missing OR multiallelic OR multi-ALT values
        ## 1 = homozygous reference (0/0) 
        ## 2 = heterozygous (0/1 or 1/0) 
        ## 3 = homozygous alternate (risk) allele (1/1)

        diploid <- lapply(strsplit(callMatrix[,1], ""), length) == 3 
        if (any(diploid == FALSE)) {
            warning("only diploid calls are supported, others are set to NA")
            callMatrix[!diploid] <- ".|." 
        }

        reflen <- width(values(mapGRanges)[["REF"]]) == 1
        altlen <- nchar(strsplit(unlist(values(mapGRanges)[["ALT"]]), ",")) == 1 
        snv <- reflen & altlen 
        if (any(snv == FALSE)) {
            warning("variants other than single nucleotide variations are ",
                    "set to NA")
            callMatrix[!snv] <- ".|." 
        }

        altelt <- elementLengths(strsplit(unlist(values(mapGRanges)[["ALT"]]), ",")) == 1 
        if (any(altelt == FALSE)) {
            warning("variants with >1 ALT allele are set to NA")
            callMatrix[!altelt] <- ".|." 
        }

        nsnps <- nrow(callMatrix)
        nsamp <- ncol(callMatrix)
        mat <-  matrix(as.raw(0), nrow=nsamp, ncol=nsnps)
        str <- strsplit(callMatrix, "")
        alt <- do.call(rbind, str)[,c(1,3)]
        alt[alt == "."] <- NA
        geno <- 2 - as.numeric(alt[,1] == 0) + as.numeric(alt[,2] == 0)
        geno[is.na(geno)] <- -1
        glst <- split(as.raw(geno+1), rep(seq(1:nsnps), nsamp))
        for (i in 1:nsnps) 
            mat[,i] = glst[[i]]

        dimnames(mat) <- list(rownames=colnames(callMatrix),
                              colnames=rownames(callMatrix))
        genotypes <- new("SnpMatrix", mat)
        flt <- !diploid | !snv | !altelt
        map <- .createMap(mapGRanges, flt)

        list(genotypes = genotypes, map = map)
    }
)

.createMap <- function(x, flt)
{
    if (is.null(x))
        DataFrame(snp.names=character(0), 
                  allele.1=DNAStringSet(), 
                  allele.2=DNAStringSet(),
                  setNA=logical())
    else 
        DataFrame(snp.names=names(x), 
                  allele.1=values(x)[["REF"]], 
                  allele.2=values(x)[["ALT"]],
                  setNA=flt)
}


