### =========================================================================
### MatrixToSnpMatrix methods 
### =========================================================================

## Coding for snpMatrix :
## 0 = missing OR multiallelic OR multi-ALT values
## 1 = homozygous reference (0|0) 
## 2 = heterozygous (0|1 or 1|0) 
## 3 = homozygous alternate (risk) allele (1|1)

setMethod("MatrixToSnpMatrix", c("matrix", "DNAStringSet", "DNAStringSetList"),
    function(callMatrix, ref, alt, ...)
{
    .Deprecated("genotypeToSnpMatrix")
  
    ok <- suppressWarnings(require("snpStats", quietly=TRUE, 
                                   character.only=TRUE))
    ok || stop("'snpStats' required; try biocLite('snpStats')", call.=FALSE) 

    map <- setNames(sapply(c(0, 1, 2, 2, 3), as.raw),
                    c(".|.", "0|0", "0|1", "1|0", "1|1"))
    diploid <- callMatrix %in% names(map)
    if (!all(diploid)) {
        warning("non-diploid variants are set to NA")
        callMatrix[!diploid] <- ".|."
    }

    altelt <- elementLengths(alt) == 1L
    if (!all(altelt)) {
        warning("variants with >1 ALT allele are set to NA")
        callMatrix[!altelt] <- ".|."
    }

    altseq <- logical(length(alt))
    idx <- rep(altelt, elementLengths(alt))
    altseq[altelt] = width(unlist(alt))[idx] == 1L
    snv <- altseq & (width(ref) == 1L)
    if (!all(snv)) {
        warning("non-single nucleotide variations are set to NA")
        callMatrix[!snv] <- ".|."
    }

    mat <- matrix(map[callMatrix], nrow=ncol(callMatrix), ncol=nrow(callMatrix),
                  byrow=TRUE, dimnames=rev(dimnames(callMatrix)))
    genotypes <- new("SnpMatrix", mat)

    flt <- !(diploid & snv & altelt)
    map <- .createMap(rownames(callMatrix), ref, alt, flt)

    list(genotypes = genotypes, map = map)
})

.createMap <- function(nms, ref, alt, flt)
{
    if (is.null(ref))
        DataFrame(snp.names=character(0), 
                  allele.1=DNAStringSet(), 
                  allele.2=DNAStringSetList(),
                  ignore=logical())
    else 
        DataFrame(snp.names=nms, 
                  allele.1=ref, 
                  allele.2=alt,
                  ignore=flt)
}
