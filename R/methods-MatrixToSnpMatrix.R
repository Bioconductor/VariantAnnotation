### =========================================================================
### MatrixToSnpMatrix methods 
### =========================================================================

## Coding for snpMatrix :
## 0 = missing OR multiallelic OR multi-ALT values
## 1 = homozygous reference (0/0) 
## 2 = heterozygous (0/1 or 1/0) 
## 3 = homozygous alternate (risk) allele (1/1)

setMethod("MatrixToSnpMatrix", c("matrix", "DNAStringSet", "DNAStringSetList"),
    function(callMatrix, ref, alt, ...)
    {
        diploid <- lapply(strsplit(callMatrix[,1], ""), length) == 3 
        if (any(diploid == FALSE)) {
            warning("only diploid calls are supported, others are set to NA")
            callMatrix[!diploid] <- ".|." 
        }

        altelt <- elementLengths(alt) == 1
        if (any(altelt == FALSE)) {
            warning("variants with >1 ALT allele are set to NA")
            callMatrix[!altelt] <- ".|." 
        }

        altseq <- lapply(alt, function(x) sum(width(x))) == 1
        altlen <- altelt & altseq 
        reflen <- width(ref) == 1
        snv <- reflen & altlen 
        if (any(snv == FALSE)) {
            warning("variants other than single nucleotide variations are ",
                    "set to NA")
            callMatrix[!snv] <- ".|." 
        }

        nsnps <- nrow(callMatrix)
        nsamp <- ncol(callMatrix)
        mat <-  matrix(as.raw(0), nrow=nsamp, ncol=nsnps)
        str <- strsplit(callMatrix, "")
        vlu <- do.call(rbind, str)[,c(1,3)]
        vlu[vlu == "."] <- NA_character_
        geno <- rowSums(matrix(as.integer(vlu), nrow=nrow(vlu))) + 1
        #geno <- 2 - as.numeric(alts[,1] == 0) + as.numeric(alts[,2] == 0)
        geno[is.na(geno)] <- 0 
        glst <- split(as.raw(geno), rep(seq(1:nsnps), nsamp))
        for (i in 1:nsnps) 
            mat[,i] = glst[[i]]

        dimnames(mat) <- list(rownames=colnames(callMatrix),
                              colnames=rownames(callMatrix))
        genotypes <- new("SnpMatrix", mat)
        flt <- !diploid | !snv | !altelt
        map <- .createMap(dimnames(callMatrix)[1], ref, alt, flt)

        list(genotypes = genotypes, map = map)
    }
)

.createMap <- function(nms, ref, alt, flt)
{
    if (is.null(ref))
        DataFrame(snp.names=character(0), 
                  allele.1=DNAStringSet(), 
                  allele.2=DNAStringSet(),
                  ignore=logical())
    else 
        DataFrame(snp.names=nms, 
                  allele.1=ref, 
                  allele.2=alt,
                  ignore=flt)
}


