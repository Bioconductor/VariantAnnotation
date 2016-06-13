### =========================================================================
### import methods 
### =========================================================================

setMethod("import", "VcfFile", function(con, format, text, ...) {
    if (!missing(format))
        checkArgFormat(con, format)
    readVcf(con, ...)
})
