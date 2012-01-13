setMethod(scanVcfHeader, "character",
    function(file, ...) 
{
    scanBcfHeader(file, ...)
})

setMethod(scanVcfHeader, "TabixFile",
    function(file, ...)
{
    scanVcfHeader(path(file), ...)
})
