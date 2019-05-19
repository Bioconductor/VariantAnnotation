### =========================================================================
### import methods 
### =========================================================================

setMethod("import", "VcfFile", function(con, format, text, ...) {
    if (!missing(format))
        rtracklayer:::checkArgFormat(con, format)
    readVcf(con, ...)
})
