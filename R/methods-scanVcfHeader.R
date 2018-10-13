setMethod(scanVcfHeader, "missing",
    function(file, ...)
{
    VCFHeader()
})

setMethod(scanVcfHeader, "character",
    function(file, ...) 
{
    if (length(file)) {
        hdr <- scanBcfHeader(file[[1]], ...)[[1]]
        VCFHeader(hdr$Reference, hdr$Sample, hdr$Header)
    } else {
        VCFHeader()
    }
})

setMethod(scanVcfHeader, "TabixFile",
    function(file, ...)
{
    scanVcfHeader(path(file), ...)
})
