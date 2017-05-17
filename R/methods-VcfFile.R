### =========================================================================
### VcfFile, VcfFileList class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

VcfFile <-
    function(file, index, ...,
             yieldSize=NA_integer_)
{
    if (is(file, "VcfFile"))
        return(file)
    tryCatch({
        Rsamtools:::.io_check_exists(file)
    }, error=function(err) {
        stop(sprintf("VcfFile: %s", conditionMessage(err)), call.=FALSE)
    })
    if (missing(index)){
        index=paste(file, "tbi", sep=".")
        if (!file.exists(index))
            index = NA
    }
    Rsamtools:::.RsamtoolsFile(.VcfFile, file, index, yieldSize=yieldSize, ...)
}

VcfFileList <-
    function(..., yieldSize=NA_integer_)
{
    Rsamtools:::.RsamtoolsFileList(
        ..., yieldSize=yieldSize, class="VcfFile",
        classDef = VariantAnnotation::VcfFile
    )
}

### Methods
###

setGeneric("indexVcf",
           function(x, ...) standardGeneric("indexVcf"), signature="x")

setMethod("indexVcf", "character",
          function(x, ...)
{
    stopifnot(!is.na(x), length(x) > 0)
    
    if(length(x) == 1L){
        index <- Rsamtools::indexTabix(x, format="vcf4", ...)
        VcfFile(x, index)        
    }else{
        index <- sapply(x, FUN=Rsamtools::indexTabix, format="vcf4", ...)
        VcfFileList(x, index)
    }    
})

setMethod("indexVcf", "VcfFile",
          function(x, ...)
{
    index <- index(x)
    newDx = NA
    if (is.na(index) || length(index) == 0L){
        index <- paste(path(x), "tbi", sep=".")
        newDx <- index
    }
    if (!file.exists(index))
        newDx <- Rsamtools::indexTabix(path(x), format="vcf4", ...)

    if (!is.na(newDx))
        index(x) <- newDx
    x    
})

setMethod("indexVcf", "VcfFileList",
          function(x, ... )
{
    ind = index(x)
    dontExist = which(!file.exists(ind))
    if (length(dontExist) != 0L){
        fn <- function(dx, obj){
            Rsamtools::indexTabix(path(obj)[dx], format="vcf4", ...)
        }
        newDx = sapply(dontExist, FUN=fn, obj=x)
        ind[dontExist] <- newDx
        index(x) <- ind
    }
    x    
})

setMethod("seqinfo", "VcfFile",
    function(x)
{
    seqinfo(scanVcfHeader(path(x)))
})

setMethod("seqinfo", "VcfFileList",
    function(x)
{
    si = lapply(x, function(elt) as(seqinfo(elt), "GRanges"))
    grl <- as(si, "GRangesList")
    seqinfo(grl)
})
