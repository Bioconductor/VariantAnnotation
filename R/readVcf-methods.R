setMethod("readVcf", c(file="character", param="missing"),
    function(file, ..., param, raw=FALSE)
{
    vcf <- scanBcf(file, character(0))
    sampleID <- scanBcfHeader(file)[[1]]$Sample
    .VcfToSummarizedExperiment(vcf, sampleID, raw=raw)
})

setMethod("readVcf", c(file="character", param="ANY"),
    function(file, index=paste(file, "tbi", sep="."), ..., param, raw=FALSE)
{
    tf <- TabixFile(file, index=paste(file, "tbi", sep="."))
    callGeneric(tf, param=param, raw=raw)
})

setMethod("readVcf", c(file="TabixFile", param="GRanges"),
    function(file, ..., param, raw=FALSE)
{
    .readVcf(file, param=param, raw=raw)
})

setMethod("readVcf", c(file="TabixFile", param="RangedData"),
    function(file, ..., param, raw=FALSE)
{
    .readVcf(file, param=param, raw=raw)
})

setMethod("readVcf", c(file="TabixFile", param="RangesList"),
    function(file, ..., param, raw=FALSE)
{
    .readVcf(file, param=param, raw=raw)
})


.readVcf <- function(file, param, raw, ...)
{
    tbx <- scanTabix(file, param=param)
    vcf <- .parseTabix(tbx, param=param)
    sampleID <- scanBcfHeader(path(file))[[1]]$Sample
    .VcfToSummarizedExperiment(vcf, sampleID, raw=raw)
}


MatrixToSnpMatrix <- function(callMatrix, ...)
{
    nsnps <- nrow(callMatrix)
    nsamples <- ncol(callMatrix)
    mat = matrix(as.raw(0), nrow=nsamples, ncol=nsnps)

    ## 0 = missing
    ## 1 = homozyg ref
    ## 2 = heterozyg ref
    ## 3 = homozyg no ref
    nalt = strsplit(callMatrix, "")
    alts <- do.call(rbind, nalt)[, c(1,3)]
    alts[alts == "."] <- NA
    nalt <- 2 - ((alts[,1] == 0) + (alts[,2] == 0))
    nalt[is.na(nalt)] <- -1
    nalt = as.raw(nalt+1)
    naltList <- split(nalt, rep(seq(1:nsnps), nsamples))
    for (i in 1:nsnps) mat[,i] = naltList[[i]]
    rownames(mat) = colnames(callMatrix)
    colnames(mat) = rownames(callMatrix)
    new("SnpMatrix", mat)
}


