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
