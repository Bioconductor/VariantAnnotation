setMethod(scanVcfHeader, "character",
    function(file, ...) 
{
    hdr <- scanBcfHeader(file, ...)[[1]]
    VCFHeader(reference=hdr$Reference, samples=hdr$Sample, 
              header=hdr$Header) 
})

setMethod(scanVcfHeader, "TabixFile",
    function(file, ...)
{
    scanVcfHeader(path(file), ...)
})
