### =========================================================================
### VCFHeader class methods 
### =========================================================================


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor 
###

VCFHeader <-
    function(reference=character(), samples=character(), 
             header=DataFrameList(), ...)
{
    new("VCFHeader", reference=reference, samples=samples, header=header, ...)
}
 
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

setMethod("reference", "VCFHeader", 
    function(x) 
{
    slot(x, "reference")
})

setMethod("samples", "VCFHeader", 
    function(object) 
{
    slot(object, "samples")
})

setMethod("header", "VCFHeader", 
    function(x) 
{
    slot(x, "header")
})

## meta
setMethod("meta", "VCFHeader", 
    function(x) 
{
    dat <- slot(x, "header")
    nms <- c("INFO", "FORMAT", "QUAL", "FILTER", "ALT", "REF")
    dat[!names(dat) %in% nms] 
})

setReplaceMethod("meta", c("VCFHeader", "DataFrame"), 
    function(x, value)  meta(x) <- as(value, "DataFrameList")
)

setReplaceMethod("meta", c("VCFHeader", "DataFrameList"), 
    function(x, value) 
{
    dat <- slot(x, "header")
    slot(x, "header") <- c(dat[!names(dat) %in% names(value)], value) 
    validObject(x)
    x
})

## fixed (QUAL, FILTER, ALT, REF)
setMethod("fixed", "VCFHeader", 
    function(x) 
{
    fixed <- c("QUAL", "FILTER", "ALT", "REF") 
    lst <- slot(x, "header")
    lst[names(lst) %in% fixed]
})

setReplaceMethod("fixed", c("VCFHeader", "DataFrameList"), 
    function(x, value) 
{
    if (!all(names(value) %in% c("QUAL", "FILTER", "ALT", "REF")))
        stop("names for 'fixed' can be 'QUAL', 'FILTER', 'ALT' or 'REF'")
    dat <- slot(x, "header")
    slot(x, "header") <- c(dat[!names(dat) %in% names(value)], value) 
    x
})

## info 
setMethod("info", "VCFHeader", 
    function(x) 
{
    info <- slot(x, "header")$INFO
    if (is.null(info))
        info <- DataFrame()
    info
})

setReplaceMethod("info", c("VCFHeader", "DataFrame"), 
    function(x, value) 
{
    slot(x, "header")$INFO <- value
    validObject(x)
    x
})

## geno
setMethod("geno", "VCFHeader",
    function(x)
{
    geno <- slot(x, "header")$FORMAT
    if (is.null(geno))
        geno <- DataFrame()
    geno 
})

setReplaceMethod("geno", c("VCFHeader", "missing", "DataFrame"),
    function(x, i, ..., value)
{
    slot(x, "header")$FORMAT <- value
    validObject(x)
    x
})

setMethod("seqinfo", "VCFHeader",
function(x)
{
  contig <- slot(x, "header")$contig
  if (is.null(contig))
    Seqinfo()
  else Seqinfo(rownames(contig),
               seqlengths =
               if (is.null(contig$length)) NA else as.integer(contig$length),
               genome = if (is.null(contig$assembly)) NA else contig$assembly)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod(show, "VCFHeader",
    function(object)
{
    selectSome <- BiocGenerics:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("class:", class(object), "\n")

    samples <- samples(object) 
    scat("samples(%d): %s\n", samples)

    meta <- names(meta(object)) 
    scat("meta(%d): %s\n", meta)

    fixed <- names(fixed(object)) 
    scat("fixed(%d): %s\n", fixed)

    info <- rownames(info(object)) 
    scat("info(%d): %s\n", info)

    geno <- rownames(geno(object))
    scat("geno(%d): %s\n", geno)
})
