### =========================================================================
### VcfFile, VcfFileList class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

VcfFile <-
    function(file, index=paste(file, "tbi", sep="."), ...,
             yieldSize=NA_integer_)
{
    tryCatch({
        Rsamtools:::.io_check_exists(c(file, index))
    }, error=function(err) {
        stop(sprintf("VcfFile: %s", conditionMessage(err)), call.=FALSE)
    })
    Rsamtools:::.RsamtoolsFile(.VcfFile, file, index, yieldSize=yieldSize, ...)
}

VcfFileList <-
    function(..., yieldSize=NA_integer_)
{
    Rsamtools:::.RsamtoolsFileList(..., yieldSize=yieldSize, class="VcfFile")
}
